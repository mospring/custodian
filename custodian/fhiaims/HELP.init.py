# Init methods for Job object from other calculators

class NwchemJob(Job):
    """
    A basic Nwchem job.
    """

    def __init__(self, nwchem_cmd, input_file="mol.nw",
                 output_file="mol.nwout", gzipped=False,
                 backup=True, settings_override=None):
        """
        Initializes a basic NwChem job.

        Args:
            nwchem_cmd ([str]): Command to run Nwchem as a list of args. For
                example, ["nwchem"].
            input_file (str): Input file to run. Defaults to "mol.nw".
            output_file (str): Name of file to direct standard out to.
                Defaults to "mol.nwout".
            backup (bool): Whether to backup the initial input files. If True,
                the input files will be copied with a ".orig" appended.
                Defaults to True.
            gzipped (bool): Deprecated. Please use the Custodian class's
                gzipped_output option instead.
            settings_override ([dict]):
                An ansible style list of dict to override changes.
                #TODO: Not implemented yet.
        """
        self.nwchem_cmd = nwchem_cmd
        self.input_file = input_file
        self.output_file = output_file
        self.backup = backup
        self.gzipped = gzipped
        self.settings_override = settings_override


class VaspJob(Job):
    """
    A basic vasp job. Just runs whatever is in the directory. But conceivably
    can be a complex processing of inputs etc. with initialization.
    """

    def __init__(self, vasp_cmd, output_file="vasp.out",
                 stderr_file="std_err.txt", suffix="", final=True,
                 backup=True, auto_npar=True, auto_gamma=True,
                 settings_override=None, gamma_vasp_cmd=None,
                 copy_magmom=False, auto_continue=False):
        """
        This constructor is necessarily complex due to the need for
        flexibility. For standard kinds of runs, it's often better to use one
        of the static constructors. The defaults are usually fine too.

        Args:
            vasp_cmd (str): Command to run vasp as a list of args. For example,
                if you are using mpirun, it can be something like
                ["mpirun", "pvasp.5.2.11"]
            output_file (str): Name of file to direct standard out to.
                Defaults to "vasp.out".
            stderr_file (str): Name of file to direct standard error to.
                Defaults to "std_err.txt".
            suffix (str): A suffix to be appended to the final output. E.g.,
                to rename all VASP output from say vasp.out to
                vasp.out.relax1, provide ".relax1" as the suffix.
            final (bool): Indicating whether this is the final vasp job in a
                series. Defaults to True.
            backup (bool): Whether to backup the initial input files. If True,
                the INCAR, KPOINTS, POSCAR and POTCAR will be copied with a
                ".orig" appended. Defaults to True.
            auto_npar (bool): Whether to automatically tune NPAR to be sqrt(
                number of cores) as recommended by VASP for DFT calculations.
                Generally, this results in significant speedups. Defaults to
                True. Set to False for HF, GW and RPA calculations.
            auto_gamma (bool): Whether to automatically check if run is a
                Gamma 1x1x1 run, and whether a Gamma optimized version of
                VASP exists with ".gamma" appended to the name of the VASP
                executable (typical setup in many systems). If so, run the
                gamma optimized version of VASP instead of regular VASP. You
                can also specify the gamma vasp command using the
                gamma_vasp_cmd argument if the command is named differently.
            settings_override ([dict]): An ansible style list of dict to
                override changes. For example, to set ISTART=1 for subsequent
                runs and to copy the CONTCAR to the POSCAR, you will provide::

                    [{"dict": "INCAR", "action": {"_set": {"ISTART": 1}}},
                     {"file": "CONTCAR",
                      "action": {"_file_copy": {"dest": "POSCAR"}}}]
            gamma_vasp_cmd (str): Command for gamma vasp version when
                auto_gamma is True. Should follow the list style of
                subprocess. Defaults to None, which means ".gamma" is added
                to the last argument of the standard vasp_cmd.
            copy_magmom (bool): Whether to copy the final magmom from the
                OUTCAR to the next INCAR. Useful for multi-relaxation runs
                where the CHGCAR and WAVECAR are sometimes deleted (due to
                changes in fft grid, etc.). Only applies to non-final runs.
            auto_continue (bool): Whether to automatically continue a run
                if a STOPCAR is present. This is very usefull if using the
                wall-time handler which will write a read-only STOPCAR to
                prevent VASP from deleting it once it finishes
        """
        self.vasp_cmd = vasp_cmd
        self.output_file = output_file
        self.stderr_file = stderr_file
        self.final = final
        self.backup = backup
        self.suffix = suffix
        self.settings_override = settings_override
        self.auto_npar = auto_npar
        self.auto_gamma = auto_gamma
        self.gamma_vasp_cmd = gamma_vasp_cmd
        self.copy_magmom = copy_magmom
        self.auto_continue = auto_continue


class QchemJob(Job):
    """
    A basic QChem Job.
    """

    def __init__(self, qchem_cmd, input_file="mol.qcinp",
                 output_file="mol.qcout", chk_file=None, qclog_file=None,
                 gzipped=False, backup=True, alt_cmd=None,
                 large_static_mem=False, total_physical_memory=100):
        """
        This constructor is necessarily complex due to the need for
        flexibility. For standard kinds of runs, it's often better to use one
        of the static constructors.

        Args:
            qchem_cmd ([str]): Command to run QChem as a list args (without
                input/output file name). For example: ["qchem", "-np", "24"]
            input_file (str): Name of the QChem input file.
            output_file (str): Name of the QChem output file.
            chk_file (str): Name of the QChem check point file. None means no
                checkpoint point file. Defaults to None.
            qclog_file (str): Name of the file to redirect the standard output
                to. None means not to record the standard output. Defaults to
                None.
            gzipped (bool): Whether to gzip the final output. Defaults to False.
            backup (bool): Whether to backup the initial input files. If True,
                the input files will be copied with a ".orig" appended.
                Defaults to True.
            alt_cmd (dict of list): Alternate commands.
                For example: {"openmp": ["qchem", "-seq", "-nt", "24"]
                              "half_cpus": ["qchem", "-np", "12"]}
            large_static_mem: use ultra large static memory
            total_physical_memory (int): The total physical memory available to
                the QChem job in unit of GB.
        """
        self.qchem_cmd = self._modify_qchem_according_to_version(copy.deepcopy(qchem_cmd))
        self.input_file = input_file
        self.output_file = output_file
        self.chk_file = chk_file
        self.qclog_file = qclog_file
        self.gzipped = gzipped
        self.backup = backup
        self.current_command = self.qchem_cmd
        self.current_command_name = "general"
        self.total_physical_memory = total_physical_memory
        self.large_static_mem = large_static_mem
        self.alt_cmd = {k: self._modify_qchem_according_to_version(c)
                        for k, c in copy.deepcopy(alt_cmd).items()}
        available_commands = ["general"]
        if self.alt_cmd:
            available_commands.extend(self.alt_cmd.keys())
        self._set_qchem_memory()


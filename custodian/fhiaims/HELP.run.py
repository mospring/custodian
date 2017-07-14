# run methods for Job objects from other codes

# NwchemJob:
    def run(self):
        """
        Performs actual nwchem run.
        """
        with zopen(self.output_file, 'w') as fout:
            return subprocess.Popen(self.nwchem_cmd + [self.input_file],
                                    stdout=fout)

# Vasp:
    def run(self):
        """
        Perform the actual VASP run.

        Returns:
            (subprocess.Popen) Used for monitoring.
        """
        cmd = list(self.vasp_cmd)
        if self.auto_gamma:
            vi = VaspInput.from_directory(".")
            kpts = vi["KPOINTS"]
            if kpts.style == Kpoints.supported_modes.Gamma \
                    and tuple(kpts.kpts[0]) == (1, 1, 1):
                if self.gamma_vasp_cmd is not None and which(
                        self.gamma_vasp_cmd[-1]):
                    cmd = self.gamma_vasp_cmd
                elif which(cmd[-1] + ".gamma"):
                    cmd[-1] += ".gamma"
        logger.info("Running {}".format(" ".join(cmd)))
        with open(self.output_file, 'w') as f_std, \
                open(self.stderr_file, "w", buffering=1) as f_err:
            # use line buffering for stderr
            p = subprocess.Popen(cmd, stdout=f_std, stderr=f_err)
        return p

# Qchem:
    def run(self):
        if "NERSC_HOST" in os.environ and (os.environ["NERSC_HOST"] in ["cori", "edison"]):
            nodelist = os.environ["QCNODE"]
            num_nodes = len(nodelist.split(","))
            tmp_creation_cmd = shlex.split("srun -N {} --ntasks-per-node 1 --nodelist {}  mkdir /dev/shm/eg_qchem".format(num_nodes, nodelist))
            tmp_clean_cmd = shlex.split("srun -N {} --ntasks-per-node 1 --nodelist {} rm -rf /dev/shm/eg_qchem".format(num_nodes, nodelist))
        elif "NERSC_HOST" in os.environ and os.environ["NERSC_HOST"] == "matgen":
            nodelist = os.environ["QCNODE"]
            num_nodes = len(nodelist.split(","))
            tmp_creation_cmd = shlex.split("mpirun -np {} --npernode 1 --host {}  mkdir /dev/shm/eg_qchem".format(num_nodes, nodelist))
            tmp_clean_cmd = shlex.split("mpirun -np {} --npernode 1 --host {} rm -rf /dev/shm/eg_qchem".format(num_nodes, nodelist))
        else:
            tmp_clean_cmd = None
            tmp_creation_cmd = None
        logging.info("Scratch dir creation command is {}".format(tmp_creation_cmd))
        logging.info("Scratch dir deleting command is {}".format(tmp_clean_cmd))
        if self.qclog_file:
            with open(self.qclog_file, "a") as filelog:
                if tmp_clean_cmd:
                    filelog.write("delete scratch before running qchem using command {}\n".format(tmp_clean_cmd))
                    subprocess.call(tmp_clean_cmd, stdout=filelog)
                if tmp_creation_cmd:
                    filelog.write("Create scratch dir before running qchem using command {}\n".format(tmp_creation_cmd))
                    subprocess.call(tmp_creation_cmd, stdout=filelog)
                returncode = self._run_qchem(log_file_object=filelog)
                if tmp_clean_cmd:
                    filelog.write("clean scratch after running qchem using command {}\n".format(tmp_clean_cmd))
                    subprocess.call(tmp_clean_cmd, stdout=filelog)
        else:
            if tmp_clean_cmd:
                subprocess.call(tmp_clean_cmd)
            if tmp_creation_cmd:
                subprocess.call(tmp_creation_cmd)
            returncode = self._run_qchem()
            if tmp_clean_cmd:
                subprocess.call(tmp_clean_cmd)
        return returncode


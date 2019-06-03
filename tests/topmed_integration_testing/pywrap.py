import sys
import subprocess


class TopMedTest():
    def __init__(self, jobname, cmd):
        self.jobname = jobname
        self.execute_program(cmd)

    def execute_program(self, cmd: str):
        """
        Execute the provided command line and check the output.
        :param cmd: str, the command to run
        """
        with open(f'{self.jobname}.stderr.txt', 'wb') as stderr_fh, \
                open(f'{self.jobname}.stdout.txt', 'wb') as stdout_fh, \
                open(f'{self.jobname}.log', 'wb') as logged:
            logged.write((cmd+'\n\n').encode('utf-8'))
            try:
                p = subprocess.run(cmd, shell=True, check=True,
                                   stdout=stdout_fh, stderr=stderr_fh, universal_newlines=True)
                logged.write((str(p)+'\n').encode('utf-8'))
            except subprocess.CalledProcessError as err:
                stderr_fh.flush()
                stdout_fh.flush()
                msg = f"\nERROR: {str(err)}\n".encode('utf-8')
                stderr_fh.write(msg)
                stdout_fh.write(msg)
                logged.write((str(err)+'\n').encode('utf-8'))
                raise err


if __name__ == "__main__":
    TopMedTest(sys.argv[1], sys.argv[2])

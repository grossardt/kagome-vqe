import signal
import time
from datetime import datetime
from qiskit_ibm_runtime import QiskitRuntimeService, Session, Estimator
from qiskit_aer.primitives import Estimator as AerEstimator
from qiskit.providers import JobStatus
from qiskit.providers.fake_provider import FakeGuadalupeV2
from qiskit_aer.noise import NoiseModel

"""Estimator classes for VQE."""

def guadalupe_noise(noise=True):
    """Returns noise model, coupling map, and basis gates for Guadalupe."""
    fake_backend = FakeGuadalupeV2()
    noise_model  = NoiseModel.from_backend(fake_backend)
    coupling_map = fake_backend.coupling_map
    basis_gates  = fake_backend.operation_names
    return {
        'coupling_map': (None,coupling_map)[noise],
        'noise_model': (None,noise_model)[noise],
        'basis_gates': (None,basis_gates)[noise],
    }

def load_qiskit_runtime_service():
    """Loads the runtime service, returns tuple (service, backend)"""
    service = QiskitRuntimeService(channel='ibm_quantum')
    backend = service.backends(simulator=True)[0]
    print('Selected backend is:', backend.name)
    return (service, backend)

def timeout_handler(signum, frame):
    raise Exception('Iteration timed out')
    
class RetryEstimator(Estimator):
    """RuntimeRetryEstimator class.
    
    This class inherits from Qiskit IBM Runtime's Estimator and overwrites its run method such that
    it retries calling it a maximum of 'max_retries' consecutive times, if it encounters one of the
    following randomly occuring errors:
    
    * An Estimator error (in this case "Job.ERROR" is printed, and the job is cancelled
      automatically)
    * A timeout error where the job either remains running or completes but does not return
      anything, for a time larger than 'timeout' (in this case the job is cancelled by the patch
      and "Job.CANCELLED" is printed)
    * A creation error, where the job fails to be created because connection is lost between the
      runtime server and the quantum computer (in this case "Failed to create job." is printed).
      If this error occurs, the patch connects the user to a new Session (to be handled with care!
      also, this will unfortunately put the next job in the queue). 
    """
    
    def __init__(
        self, *args, noise: bool = False, max_retries: int = 5, timeout: int = 3600, **kwargs
    ) -> None:
        super().__init__(*args, **kwargs)
        self.set_options(
            simulator=guadalupe_noise(noise), 
            transpilation={'skip_transpilation': True},
        )
        self.max_retries = max_retries
        self.timeout = timeout
        self.backend = super().session._backend
        signal.signal(signal.SIGALRM, timeout_handler)
    
    def run(self, circuits, observables, parameter_values, **kwargs):
        result = None
        for i in range(self.max_retries):
            try:
                job = super().run(circuits, observables, parameter_values, **kwargs)
                while job.status() in [JobStatus.INITIALIZING, JobStatus.QUEUED,
                                       JobStatus.VALIDATING]:
                    print(datetime.now().strftime('%Y%m%d-%H%M%S: '), end='')
                    print('Estimator attempt {} of {}, status: {}'.format(
                        i+1, self.max_retries, job.status().value), end='\r')
                    time.sleep(5) # Check every 5 seconds whether job status has changed
                signal.alarm(self.timeout) # Once job starts running, set timeout to 1h by default
                result = job.result()
                if result is not None:
                    signal.alarm(0) # reset timer
                    return job
            except Exception as e:
                print("\nSomething went wrong...")
                print(f"\n\nERROR MESSAGE:\n{e}\n\n")
                if 'job' in locals(): # Sometimes job fails to create
                    print(f"Job ID: {job.job_id}. Job status: {job.status()}.")
                    if job.status() not in [JobStatus.DONE, JobStatus.ERROR, JobStatus.CANCELLED]:
                        job.cancel()
                else:
                    print("Failed to create job.")
                print(f"Starting trial number {i+2}...\n")
                print(f"Creating new session...\n")
                signal.alarm(0) # reset timer
                super().session.close()
                self._session = Session(backend=self.backend)
        if result is None:
            raise RuntimeError(
                f"Program failed! Maximum number of retries ({self.max_retries}) exceeded")
        
class LocalEstimator(AerEstimator):
    """Estimator for local simulation.
    
    This class inherits from Qiskit's AerAestimator and overwrites its constructor such that it can
    take a simplified list of parameters:
    LocalEstimator(..., noise = True/False, gpu = True/False, ...)
    sets the correct backend options for the AerSimulator to run with or without noise
    and on CPU or GPU, respectively.
    It can take all other parameters of the AerEstimator constructor.
    """
    
    def __init__(self, *args, noise: bool = False, gpu: bool = False, **kwargs):
        backend_options = {
                'method': 'density_matrix',
                'device': ('CPU','GPU')[gpu],
            }
        backend_options.update(guadalupe_noise(noise))
        super().__init__(
            *args,
            backend_options=backend_options,
            run_options={'shots': 1024},
            skip_transpilation=True,
            **kwargs,
        )

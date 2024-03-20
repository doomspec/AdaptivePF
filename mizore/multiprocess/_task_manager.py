from multiprocessing import Process, Queue
import time
import numpy
from copy import deepcopy
from ._task import Task
from tqdm import tqdm

RECEIVE_PERIOD = 0.01

class TaskResult:

    def __init__(self, result, index_of_in, series_id):
        self.index_of_in = index_of_in
        self.result = result
        self.series_id = series_id


class TaskRunner(Process):
    """
    See TaskManager
    """

    def __init__(self, task_queue, result_queue):
        Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        while True:
            #time.sleep(RECEIVE_PERIOD)
            task_package = self.task_queue.get(True)
            # print(self.name,"get task")
            result_package = []
            public_resource = task_package[0]
            for i in range(1, len(task_package)):
                task = task_package[i]
                dress_by_public_resource(public_resource, task)
                result = task.run()
                result_package.append(TaskResult(result, task.index_of_in, task.series_id))
            self.result_queue.put(result_package)


def dress_by_public_resource(public_resource: dict, obj):
    if public_resource is None:
        return
    for key in public_resource.keys():
        obj.__dict__[key] = public_resource[key]


class TaskManager:
    """
    The class for parallel evaluate quantum circuits.
    Holds the task queue and result queue for TaskRunner to read the inputs and put the outputs.
    Tasks are distributed to different TaskRunner running on different cores to implement parallel
    Usage:
        1. Add tasks to the buffer and attach a *task_series_id*: add_task_to_buffer(_task, task_series_id)
        2. Use flush() to send the tasks in the buffer to the TaskRunners
        3. Use receive_task_result(task_series_id) to receive the results of the tasks with *task_series_id*
    Advanced Usage:
        When flush(), a *public_resource* can be added to avoid including 
        large and common resources (like a big Hamiltonian) in every tasks
        public_resource should be a dict() with like {"hamiltonian":operator}
        By doing so, the TaskRunner will replace the attribute named "hamiltonian" by operator in the tasks before run 
    """

    def __init__(self, n_processor=4, task_package_size=5):

        self.n_processor = n_processor
        self.processor_list = []
        self.task_queue = Queue()
        self.result_queue = Queue()
        self.n_task_processed = 0
        self.n_task_remain_by_series_id = dict()
        self.buffer_to_send = []
        self.receive_buffer_by_series_id = dict()
        self.task_package_size = task_package_size

        self.task_difficulty_dict = dict()

        for i in range(n_processor):
            self.processor_list.append(TaskRunner(
                self.task_queue, self.result_queue))
        for i in range(n_processor):
            self.processor_list[i].start()

    def add_task_to_buffer(self, _task: Task, task_series_id=0):
        task = deepcopy(_task)
        task.index_of_in = self.n_task_processed
        task.series_id = task_series_id
        self.n_task_processed += 1
        if task_series_id in self.n_task_remain_by_series_id.keys():
            self.n_task_remain_by_series_id[task_series_id] += 1
        else:
            self.n_task_remain_by_series_id[task_series_id] = 1
            self.receive_buffer_by_series_id[task_series_id] = []
        self.buffer_to_send.append(task)

    def smart_flush(self,public_resource=None,mode="average"):
        task_package_size=(len(self.buffer_to_send)//self.n_processor)+1
        self.flush(task_package_size=task_package_size,public_resource=public_resource)

    def flush(self, task_package_size=None, public_resource=None):
        if task_package_size is None:
            task_package_size = self.task_package_size

        task_package = [public_resource]
        for task in self.buffer_to_send:
            task_package.append(task) # TODO avoid using append!
            if len(task_package) >= task_package_size:
                self.task_queue.put(task_package)
                task_package = [public_resource]
        if len(task_package) != 0:
            self.task_queue.put(task_package)
        self.buffer_to_send = []


    def receive_task_result(self, task_series_id=0, progress_bar=False):

        result_list = []
        index_list = []
        n_total_tasks = self.n_task_remain_by_series_id[task_series_id]
        for task_result in self.receive_buffer_by_series_id[task_series_id]:
            result_list.append(task_result.result)
            index_list.append(task_result.index_of_in)
            self.n_task_remain_by_series_id[task_series_id] -= 1
        self.receive_buffer_by_series_id[task_series_id] = []

        if progress_bar:
            pbar = tqdm(total=n_total_tasks)
            pbar.set_description(str(task_series_id))
            last_progress_value = n_total_tasks

        while self.n_task_remain_by_series_id[task_series_id] != 0:
            result_package = self.result_queue.get(True)
            for task_result in result_package:
                if task_result.series_id == task_series_id:
                    result_list.append(task_result.result)
                    index_list.append(task_result.index_of_in)
                    self.n_task_remain_by_series_id[task_series_id] -= 1
                else:
                    self.receive_buffer_by_series_id[task_result.series_id].append(task_result)

            if progress_bar:
                pbar.update(last_progress_value - self.n_task_remain_by_series_id[task_series_id])
                last_progress_value = self.n_task_remain_by_series_id[task_series_id]

        if progress_bar:
            pbar.close()

        index_rank_list = numpy.argsort(numpy.array(index_list))
        ranked_result_list = []

        for i in range(len(index_list)):
            ranked_result_list.append(result_list[index_rank_list[i]])
        return ranked_result_list

    def close(self):
        for processor in self.processor_list:
            processor.terminate()

    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

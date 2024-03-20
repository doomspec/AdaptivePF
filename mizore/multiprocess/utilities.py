


def get_task_packages(task_list,task_package_size=10):
    n_task=len(task_list)
    n_package=n_task//task_package_size

    task_package_list=[0]*n_package
    for i in range(n_package):
        task_package=[0]*task_package_size
        for j in range(task_package_size):
            task_package[j]=task_list[i*task_package_size+j]
        task_package_list[i]=task_package

    n_task_remain=n_task%task_package_size

    if n_task_remain==0:
        return task_package_list
    else:
        task_package = [0] * n_task_remain
        for j in range(n_task_remain):
            task_package[j]=task_list[n_task-n_task_remain+j]
        task_package_list.append(task_package)
        return task_package_list
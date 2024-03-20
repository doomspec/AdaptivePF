from heapq import nlargest,nsmallest


def get_top_k_indices(scores, k):
    return nlargest(k, range(len(scores)), scores.__getitem__)

def get_bottom_k_indices(scores, k):
    return nsmallest(k, range(len(scores)), scores.__getitem__)

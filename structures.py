__author__ = 'nishaswarup'

def heap_push(heap, value):
    '''
    This function pushes the element 'value' into the heap
    By pushing it up through the heap to it's proper location.
    This is for a max-value heap.
    :param heap: The heap we're adding to. We assume it's properly formatted.
    :param value: The value we want to add to the heap
    :return: None
    '''
    # we add a temporary variable to the heap
    heap.append(0)
    value_pos = len(heap) - 1
    # this is the parent's location
    temp = max((value_pos + 1)/2 - 1, 0)
    # a node's parent is located at pos/2
    # we push value through the heap, pulling everything below it
    # down one level in the heap
    while heap[temp] < value:
        heap[value_pos] = heap[temp]
        value_pos = temp
        # we break if the elmeent is in the beginning of the array
        if value_pos == 0:
            break

        temp = max((value_pos + 1)/2 - 1, 0)

    heap[value_pos] = value

def heap_pop(heap):
    '''
    This function pops the largest value from our max-value heap,
    reconstructs the heap, and returns the former largest value
    :param heap: A properly sorted max-value heap
    :return: The max element in heap, it gets removed
    '''
    # handle the edge case if heap is of size 0
    if len(heap) == 1:
        return heap.pop()

    # we want to pop the last element,
    # because that's much faster
    out = heap[0]
    heap[0] = heap.pop()
    # now we move the last element down the list
    # to it's proper location
    down_pos = 0
    while down_pos < len(heap):
        child1 = down_pos*2 + 1
        child2 = child1 + 1
        bigger_child = child1
        if child1 >= len(heap): break
        if child2 < len(heap) and heap[child2] > heap[child1]:
            bigger_child = child2
        # now bigger child has the index for the larger of the two children
        # if the head node is smaller than its biggest child, we swap
        if heap[down_pos] < heap[bigger_child]:
            heap[bigger_child], heap[down_pos] = heap[down_pos], heap[bigger_child]
            down_pos = bigger_child
        # if the element is in the right spot, we return
        else:
            break
    return out
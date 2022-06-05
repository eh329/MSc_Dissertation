def array_change(arr):

    make = 0

    for i in range(0, len(arr) - 1):
        if arr[i + 1] > arr[i]:
            continue
        
        elif arr[i + 1] == arr[i]:
            
            if i == 0:
                make += 1

            else:
                make += 1 + i
            
        else:
            make += (abs((arr[i + 1] - arr[i])) + 1)
            arr[i + 1] = arr[i + 1] + make

    return make
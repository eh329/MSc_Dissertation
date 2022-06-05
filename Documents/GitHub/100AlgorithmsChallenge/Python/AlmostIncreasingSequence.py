def increasing_sequence(arr):
    mistakes = 0
    
    for i in range(0, len(arr) - 1):
        if arr[i] + 1 > arr[i + 1]:
            mistake += 1
        
        else:
            continue


    if mistake > 1:
        return False
    
    else:
        return True   




## 1. Drawdown Calculation

Based on returns:

    # 1. Convert returns to a cumulative wealth (price)
    wealth_index = (1 + returns).cumprod()
    
    # 2. Calculate the expanding maximum (running peak)
    running_max = wealth_index.cummax()
    
    # 3. Calculate the percentage drawdown at each point
    drawdowns = (running_max - wealth_index) / running_max
    
    # 4. Return the absolute maximum drawdown
    return drawdown.max()

Based on price:

    # 1. Calculate the expanding maximum (running peak)
    running_max = prices.cummax()
    
    # 2. Calculate the percentage drawdown at each point
    # Formula: (Peak - Current Price) / Peak
    drawdowns = (running_max - prices) / running_max
    
    # 3. Return the absolute maximum drawdown
    return drawdowns.max()

## 2. Sorting Algorithm

def quick_sort(arr):

    if len(arr) <= 1:

        return arr

    pivot = arr[0]

    left = [x for x in arr if x < pivot]

    middle = [x for x in arr if x == pivot]

    right = [x for x in arr if x > pivot]

    return quick_sort(left) + middle + quick_sort(right)

  
import pandas as pd
import numpy as np

left={
    'run_id': np.repeat([30, 15,], 10),
    'array_event_id' : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12],
    'something' : np.random.normal(size=20),
}


right={
    'run_id': np.repeat([30, 15,], 7),
    'array_event_id' : [1, 2, 3, 7, 8, 9, 10, 2, 3, 4, 5, 6, 9, 12],
    'something' : np.random.normal(size=14),
}

left = pd.DataFrame(left).set_index(['run_id', 'array_event_id'])
right = pd.DataFrame(right).set_index(['run_id', 'array_event_id'])

print(left)
print(right)

print('MErge Result')
print(pd.merge(left, right, left_index=True, right_index=True))
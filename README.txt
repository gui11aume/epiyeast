Each segment is in a different directory.
To reproduce the data, enter a directory and type "make".

Preprocessing steps.

1. The reads are oriented using the constant part.
2. The reads are demultiplexed from the primer information.
3. The frame is found by finding the constant part on the second read.
4. The reads are merged, resolving conflicts with the quality.

Reads that fail any of those steps are ignored.

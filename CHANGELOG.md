### What’s new in v1.1.0
1. Each new segment is connected to the last segment. We do this by replacing the first point of the new segment with the last point of the previous segment.
2. If the distance between the last point of the previous segment and the first point of the new segment is greater than Delta_min, then print a WARNING.

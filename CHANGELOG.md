### What’s new in v2.2.0
1. Each new segment is connected to the last segment. We do this by replacing the first point of the new segment with the last point of the previous segment.
   
2. If the distance between the last point of the previous segment and the first point of the new segment is greater than Delta_min, then print a WARNING.

3. Change in the syntax of the structure. They are now manif.points.pos instead of manif.pointspos and manif.points.neg instead of manif.pointsneg. If the manifold is orientation-preserving, then anyway we call them .pos or .neg, depending on which direction we have computed the manifold.

4. Changes in manifplot and eps_orbit, to be consistent with the new syntax.

5. Add a function to add a new branch to the structure.

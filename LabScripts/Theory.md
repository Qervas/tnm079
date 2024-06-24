[toc]

# lab 1

## 1. mean curvature

### Mean Curvature

The mean curvature $ H $ at a vertex $ v_i $ can be estimated using a discrete approximation based on the cotangent formula and the Voronoi area. The mean curvature vector $ \mathbf{H}_i $ can be given by:

$ \mathbf{H}_i = \frac{1}{4A_i} \sum_{j \in N_1(i)} (\cot \alpha_j + \cot \beta_j)(\mathbf{v}_i - \mathbf{v}_j) $

Here:

- $ A_i $ is the Voronoi area associated with vertex $ v_i $.
- $ N_1(i) $ represents the 1-ring neighborhood of $ v_i $.
- $ \alpha_j $ and $ \beta_j $ are the angles opposite to the edge $ (\mathbf{v}_i, \mathbf{v}_j) $ in the triangles sharing this edge.

### Voronoi Area

The Voronoi area $ A_i $ associated with vertex $ v_i $ is an important component for normalizing the curvature vector. It can be computed as:

$ A_i = \frac{1}{8} \sum_{j \in N_1(i)} (\cot \alpha_j + \cot \beta_j) \| \mathbf{v}_i - \mathbf{v}_j \|^2 $

### Detailed Steps

1. **Calculate Cotangent Weights**: For each vertex $ v_i $, find its neighboring vertices $ v_j $. For each pair of adjacent neighbors $ v_j $ and $ v_k $, calculate the cotangent of the angles $ \alpha_j $ and $ \beta_j $.
2. **Compute Voronoi Area**: Use the cotangent weights to calculate the Voronoi area $ A_i $ for vertex $ v_i $.
3. **Estimate Mean Curvature Vector**: Sum up the contributions from all neighbors $ v_j $ using the cotangent weights and the differences in vertex positions.
4. **Normalize Mean Curvature Vector**: Divide the summed vector by $ 4A_i $ to get the mean curvature vector $ \mathbf{H}_i $.
5. **Calculate Mean Curvature**: The scalar mean curvature $ H $ is the magnitude of the mean curvature vector $ \mathbf{H}_i $.

### Mathematical Formulation

1. **Cotangent Weights**:
   $ \cot \alpha_j = \frac{\mathbf{a} \cdot \mathbf{b}}{\| \mathbf{a} \times \mathbf{b} \|} $
   $ \cot \beta_j = \frac{\mathbf{c} \cdot \mathbf{d}}{\| \mathbf{c} \times \mathbf{d} \|} $
2. **Voronoi Area**:
   $ A_i = \frac{1}{8} \sum_{j \in N_1(i)} (\cot \alpha_j + \cot \beta_j) \| \mathbf{v}_i - \mathbf{v}_j \|^2 $
3. **Mean Curvature Vector**:
   $ \mathbf{H}_i = \frac{1}{4A_i} \sum_{j \in N_1(i)} (\cot \alpha_j + \cot \beta_j)(\mathbf{v}_i - \mathbf{v}_j) $
4. **Mean Curvature**:
   $ H_i = \| \mathbf{H}_i \| $

### Example Calculation

1. **Identify 1-Ring Neighbors**: For a vertex $ v_i $, find its 1-ring neighbors $ v_j $.
2. **Compute Cotangent Weights**:

   - For each neighbor $ v_j $ and the next neighbor $ v_k $ in the loop, compute the angles $ \alpha_j $ and $ \beta_j $ opposite to the edges.
3. **Calculate Voronoi Area**:

   - Sum up the cotangent weights and the squared distances.
4. **Compute Mean Curvature Vector**:

   - Use the cotangent weights to sum the weighted differences $ (\mathbf{v}_i - \mathbf{v}_j) $.
5. **Normalize and Compute Mean Curvature**:

   - Divide by $ 4A_i $ to get the mean curvature vector and then find its magnitude.

## 2. Genus

Euler characteristic:

$ \chi = V - E + F $

where $ V $ is the number of vertices, $ E $ is the number of edges, and $ F $ is the number of faces. For a closed manifold mesh with one shell (genus $ g $), the relationship between the Euler characteristic and the genus is given by:

$ \chi = 2 - 2g $

This calculates the genus correctly for a closed manifold mesh with one shell by returning:

$ g = \frac{V - E + F}{2} $

### Explanation of the Genus Calculation

1. **Vertices (V)**: These are the points in 3D space that define the corners of the mesh's faces.
2. **Edges (E)**: These are the lines connecting the vertices.
3. **Faces (F)**: These are the polygonal regions bounded by edges, usually triangles or quadrilaterals in a mesh.

The **Euler characteristic** $\chi$ for a mesh is calculated as:

$ \chi = V - E + F $

For a mesh that is a closed surface, the relationship between the Euler characteristic and the genus $ g $ is:

$ \chi = 2 - 2g $

Solving for $ g $:

$ g = \frac{2 - \chi}{2} $

By rearranging the Euler characteristic formula, we obtain the genus directly as:

$ g = \frac{V - E + F}{2} $

Genus represents the number of "holes" or handles in the mesh.

# lab 2

## 1. Maths formula derivation

Derivation of the expression:

1. The Quadric Error Function

The error function for a vertex $ v $ is given by:

$\Delta = v^T Q v$

where:

- $ v $ is the homogeneous coordinate vector $[x, y, z, 1]^T$
- $ Q $ is the 4x4 symmetric matrix representing the quadric

### 2. Partial Derivatives of the Error Function

To minimize $ \Delta $, we need to set the partial derivatives of $ \Delta $ with respect to $ x, y, $ and $ z $ to zero:

$\frac{\partial \Delta}{\partial x} = \frac{\partial \Delta}{\partial y} = \frac{\partial \Delta}{\partial z} = 0$

### 3. Expressing $\Delta$  in Terms of $ x, y, z $

Expanding $ $\Delta = v^T Q v $$:

$\Delta = \begin{pmatrix} x & y & z & 1 \end{pmatrix} \begin{pmatrix} q_{11} & q_{12} & q_{13} & q_{14} \\ q_{12} & q_{22} & q_{23} & q_{24} \\ q_{13} & q_{23} & q_{33} & q_{34} \\ q_{14} & q_{24} & q_{34} & q_{44} \end{pmatrix} \begin{pmatrix} x \\ y \\ z \\ 1 \end{pmatrix}$

This can be written as:

$\Delta = q_{11}x^2 + 2q_{12}xy + 2q_{13}xz + 2q_{14}x + q_{22}y^2 + 2q_{23}yz + 2q_{24}y + q_{33}z^2 + 2q_{34}z + q_{44}$

### 4. Taking Partial Derivatives

Taking the partial derivatives with respect to $ x, y, $ and $ z $:

$\frac{\partial \Delta}{\partial x} = 2q_{11}x + 2q_{12}y + 2q_{13}z + 2q_{14} = 0$
$\frac{\partial \Delta}{\partial y} = 2q_{12}x + 2q_{22}y + 2q_{23}z + 2q_{24} = 0$
$\frac{\partial \Delta}{\partial z} = 2q_{13}x + 2q_{23}y + 2q_{33}z + 2q_{34} = 0$

Dividing each equation by 2, we get:

$q_{11}x + q_{12}y + q_{13}z + q_{14} = 0$
$q_{12}x + q_{22}y + q_{23}z + q_{24} = 0$
$q_{13}x + q_{23}y + q_{33}z + q_{34} = 0$

### 5. Matrix Form

These three equations can be expressed in matrix form:

$\begin{pmatrix} q_{11} & q_{12} & q_{13} \\ q_{12} & q_{22} & q_{23} \\ q_{13} & q_{23} & q_{33} \end{pmatrix} \begin{pmatrix} x \\ y \\ z \end{pmatrix} = -\begin{pmatrix} q_{14} \\ q_{24} \\ q_{34} \end{pmatrix}$

### 6. Solving for $[x, y, z]^T$

If the 3x3 submatrix (top-left 3x3 part of $ Q $) is invertible, we can solve for $[x, y, z]^T$:

$\begin{pmatrix} x \\ y \\ z \end{pmatrix} = - \begin{pmatrix} q_{11} & q_{12} & q_{13} \\ q_{12} & q_{22} & q_{23} \\ q_{13} & q_{23} & q_{33} \end{pmatrix}^{-1} \begin{pmatrix} q_{14} \\ q_{24} \\ q_{34} \end{pmatrix}$

This gives us the optimal position for the new vertex after collapsing an edge.

## 2. Heuristic Quadratic error

**Objective:** Improve the mesh decimation process by introducing a heuristic that considers the camera's position, ensuring that areas closer to the camera are preserved with higher detail.

#### Quadric Error Metric

The quadric error metric (QEM) is a common method for mesh decimation. It calculates the error introduced by collapsing an edge into a single vertex. The error is derived from the sum of squared distances from the new vertex position to the planes of all faces originally incident to the two vertices being collapsed.

Given a plane equation:
$ ax + by + cz + d = 0 $

The corresponding quadric is:
$ K = \begin{pmatrix}
a^2 & ab & ac & ad \\
ab & b^2 & bc & bd \\
ac & bc & c^2 & cd \\
ad & bd & cd & d^2 \\
\end{pmatrix} $

For a vertex $ v $, the error is:
$ \text{Error}(v) = v^T Q v $
where $ Q $ is the sum of the quadrics of all planes incident to $ v $.

#### Incorporating Camera Distance Heuristic

To enhance this metric by considering the camera's position, we introduce a distance-based weight. The new cost function becomes a combination of the original QEM and a term that penalizes vertices based on their distance from the camera.

Let:

- $ \mathbf{p}_c $ be the camera position
- $ \mathbf{p}_v $ be the vertex position
- $ d(\mathbf{p}_v, \mathbf{p}_c) $ be the Euclidean distance between the vertex and the camera

The adjusted cost function can be written as:
$ \text{Cost}(v) = \mathbf{v}^T Q \mathbf{v} + \lambda \cdot d(\mathbf{p}_v, \mathbf{p}_c) $

Where:

- $ \lambda $ is a weight factor that controls the influence of the distance heuristic.

#### Steps for Incorporating the Heuristic

1. **Calculate Quadric Error:**

   - For vertices $ v_1 $ and $ v_2 $:
     $ \text{Error}(v_1) = \mathbf{v_1}^T Q_1 \mathbf{v_1} $
     $ \text{Error}(v_2) = \mathbf{v_2}^T Q_2 \mathbf{v_2} $
     $ \text{Error}(\text{midpoint}) = \mathbf{v}_{\text{mid}}^T (Q_1 + Q_2) \mathbf{v}_{\text{mid}} $
2. **Calculate Distance to Camera:**

   - For each vertex:
     $ d(\mathbf{p}_v, \mathbf{p}_c) = \|\mathbf{p}_v - \mathbf{p}_c\| $
3. **Combine Metrics:**

   - Adjust the cost for each vertex and midpoint:
     $ \text{Cost}(v_1) = \mathbf{v_1}^T Q_1 \mathbf{v_1} + \lambda \cdot d(\mathbf{p}_v1, \mathbf{p}_c) $
     $ \text{Cost}(v_2) = \mathbf{v_2}^T Q_2 \mathbf{v_2} + \lambda \cdot d(\mathbf{p}_v2, \mathbf{p}_c) $
     $ \text{Cost}(\text{midpoint}) = \mathbf{v}_{\text{mid}}^T (Q_1 + Q_2) \mathbf{v}_{\text{mid}} + \lambda \cdot d(\mathbf{p}_{\text{mid}}, \mathbf{p}_c) $
4. **Select Optimal Vertex Position:**

   - Choose the position with the minimal combined cost.

### Deduction

1. **Quadric Error Calculation:**

   - Each vertex's error is computed by summing the squared distances to its incident planes using the quadric $ Q $.
2. **Distance-Based Adjustment:**

   - By introducing the distance term, vertices closer to the camera incur a higher cost unless they are very significant in shape preservation, thus favoring preservation of detail where it is more visible.
3. **Heuristic Influence:**

   - The weight factor $ \lambda $ allows control over the trade-off between shape preservation and view-dependent simplification, ensuring a balanced decimation.
4. **Algorithm Efficiency:**

   - This heuristic does not significantly complicate the existing algorithm but provides a substantial improvement in visual quality by preserving detail where it matters most (i.e., near the camera).

# lab3

### 1. Error Explanation:

#### Subdivision Process

The subdivision process for uniform cubic splines generates a new set of control points by interpolating between the existing ones. The goal is to produce a refined curve that converges to the true shape of the original curve as the number of points increases.

#### Subdivision Rules

The rules for generating the new control points are derived from the properties of cubic B-splines. For a given set of control points $ \{ \mathbf{P}_i \} $, the new control points $ \{ \mathbf{Q}_j \} $ are computed as follows:

1. **New Control Points at Existing Locations**:
   $
   \mathbf{Q}_{2i} = \mathbf{P}_i
   $
   These points are directly carried over from the original set.
2. **New Control Points Between Existing Points**:
   $
   \mathbf{Q}_{2i+1} = \frac{1}{8} \left( \mathbf{P}_{i-1} + 6\mathbf{P}_i + \mathbf{P}_{i+1} \right)
   $
   This rule generates new points between each pair of existing points, using a weighted average to ensure smoothness.

#### Boundary Conditions

1. **First Boundary Point**:
   $
   \mathbf{Q}_0 = \mathbf{P}_0
   $
   The first control point is directly added to the new set of control points.

2. **Last Boundary Point**:
   $
   \mathbf{Q}_{2n-1} = \mathbf{P}_{n-1}
   $
   The last control point is directly added to the new set of control points.

3. **Second-to-Last Point**:
   $
   \mathbf{Q}_{2(n-1)+1} = \frac{1}{8} \left( \mathbf{P}_{n-2} + 6\mathbf{P}_{n-1} + \mathbf{P}_{n} \right)
   $

### Subdivision Algorithm

1. **Initialize**:
   - Create an empty vector `newc` to store the new control points.
   - Ensure there are at least 5 control points to apply the subdivision.

2. **Add First Boundary Point**:
   - Add the first control point directly to `newc`.

3. **Apply Subdivision Rules to Internal Points**:
   - Loop through the internal control points (excluding the first and last points).
   - For each internal point, calculate two new points using the subdivision rules and add them to `newc`.

4. **Handle Last Boundary Point**:
   - Add the second-to-last and the last control points to `newc` using the modified rule for the second-to-last point and directly for the last point.

5. **Verification**:
   - Ensure that the number of new control points is twice the original number minus one.



### 2.  cubic spline `GetValue`

1.**Starting Index Calculation:** We calculate the starting index of the control points by using `std::floor(t) - 1`. This ensures that we start from the control point just before the point ttt.

2.**Boundary Check:** We ensure that the starting index does not go out of the valid range of control points.

3.**Summation of Contributions:** We iterate over the four relevant control points, calculate the corresponding basis function value, and accumulate the weighted control points.

4.**Normalization:** Finally, we normalize the resulting value by dividing it by the sum of the basis function values to ensure proper interpolation

# lab4

## 1.Implementing CSG Operators

1. **CSG Boolean Operators:**
   - **Union:** Combine two implicit surfaces $A$ and $B$ using:
     $
     \text{Union}(A, B) = \min(A, B)
     $
   - **Intersection:** Combine two implicit surfaces $A$ and $B$ using:
     $
     \text{Intersection}(A, B) = \max(A, B)
     $
   - **Difference:** Subtract implicit surface $B$ from $A$ using:
     $
     \text{Difference}(A, B) = \max(A, -B)
     $

### Quadric Surface Implementation

1. **Quadric Equation:**

   The quadric surface is defined by the quadratic function:
   $
   f(x, y, z) = \mathbf{p}^T \mathbf{Q} \mathbf{p}
   $
   where:

   - $ \mathbf{p} = [x, y, z, 1]^T $
   - $ \mathbf{Q} $ is the 4x4 matrix of coefficients.
2. **Matrix Form of Quadric:**
   $
   f(x, y, z) = \begin{pmatrix} x & y & z & 1 \end{pmatrix}
   \begin{pmatrix}
   q_{11} & q_{12} & q_{13} & q_{14} \\
   q_{12} & q_{22} & q_{23} & q_{24} \\
   q_{13} & q_{23} & q_{33} & q_{34} \\
   q_{14} & q_{24} & q_{34} & q_{44}
   \end{pmatrix}
   \begin{pmatrix} x \\ y \\ z \\ 1 \end{pmatrix}
   $
3. **Expression for $ f(x, y, z) $:**
   $
   f(x, y, z) = q_{11}x^2 + 2q_{12}xy + 2q_{13}xz + 2q_{14}x + q_{22}y^2 + 2q_{23}yz + 2q_{24}y + q_{33}z^2 + 2q_{34}z + q_{44}
   $
4. **Calculating the Gradient:**

   The gradient $ \nabla f $ of the quadric surface is given by:
   $
   \nabla f(x, y, z) = 2 \mathbf{Q}_{\text{sub}} \mathbf{p}
   $
   where $ \mathbf{Q}_{\text{sub}} $ is the top-left 3x3 submatrix of $ \mathbf{Q} $:
   $
   \nabla f(x, y, z) = 2 \begin{pmatrix}
   q_{11} & q_{12} & q_{13} \\
   q_{12} & q_{22} & q_{23} \\
   q_{13} & q_{23} & q_{33}
   \end{pmatrix}
   \begin{pmatrix} x \\ y \\ z \end{pmatrix} +
   2 \begin{pmatrix} q_{14} \\ q_{24} \\ q_{34} \end{pmatrix}
   $

### Implementing the Gradient Operator for Implicits

1. **Discrete Gradient:**

   Implement the `GetGradient` function using the central difference approximation:
   $
   \frac{\partial f}{\partial x} \approx \frac{f(x + \epsilon, y, z) - f(x - \epsilon, y, z)}{2\epsilon}
   $
   $
   \frac{\partial f}{\partial y} \approx \frac{f(x, y + \epsilon, z) - f(x, y - \epsilon, z)}{2\epsilon}
   $
   $
   \frac{\partial f}{\partial z} \approx \frac{f(x, y, z + \epsilon) - f(x, y, z - \epsilon)}{2\epsilon}
   $

### Visualizing and Verifying the Implementations

1. **CSG Operations:**

   - Use the CSG panel in the GUI to apply union, intersection, and difference operations on multiple objects.
   - Adjust the "Mesh sampling" parameter to see how the resulting mesh changes.
2. **Quadric Surface Visualization:**

   - Load the quadric surface and visualize it using the gradient mode.
   - Compare visualized gradients with mesh normals to ensure accuracy.
3. **Discrete Gradient Visualization:**

   - Visualize the computed gradient vectors at the surface or on a vector cut plane.
   - Experiment with the "Differential scale" slider to see the effect of different $ \epsilon $ values on gradient computation.

## 2.**Discrete Curvature:**

- Implement the `GetCurvature` function using the finite difference approximation for second partial derivatives.

#### Gradient of an Implicit Surface

The gradient of an implicit surface $ \phi(x, y, z) $ at a point gives the direction of the steepest ascent. It is a vector that points perpendicular to the surface at that point. Mathematically, the gradient is defined as:

$$
\nabla \phi = \left( \frac{\partial \phi}{\partial x}, \frac{\partial \phi}{\partial y}, \frac{\partial \phi}{\partial z} \right)
$$

In the code, we use central differencing to approximate the partial derivatives:

$ \frac{\partial \phi}{\partial x} \approx \frac{\phi(x + \epsilon, y, z) - \phi(x - \epsilon, y, z)}{2\epsilon} $
$ \frac{\partial \phi}{\partial y} \approx \frac{\phi(y + \epsilon, y, z) - \phi(y - \epsilon, y, z)}{2\epsilon} $

$$
\frac{\partial \phi}{\partial z} \approx \frac{\phi(x, y, z + \epsilon) - \phi(x, y, z - \epsilon)}{2\epsilon}
$$

Where $ \epsilon $ is a small value (delta) used for numerical differentiation.

#### Curvature of an Implicit Surface

Curvature measures how much the surface deviates from being flat. For an implicit surface, the mean curvature can be calculated using the gradient and second-order partial derivatives. The formula for curvature $ \kappa $ is derived from differential geometry:
$ \kappa = \nabla \cdot \left( \frac{\nabla \phi}{|\nabla \phi|} \right) $

This can be expanded using the second-order partial derivatives of $ \phi $:

$$
\kappa = \frac{\phi_{xx} (\phi_y^2 + \phi_z^2) - 2 \phi_x \phi_y \phi_{xy} - 2 \phi_x \phi_z \phi_{zx} + \phi_{yy} (\phi_x^2 + \phi_z^2) - 2 \phi_y \phi_z \phi_{yz} + \phi_{zz} (\phi_x^2 + \phi_y^2)}{(\phi_x^2 + \phi_y^2 + \phi_z^2)^{3/2}}
$$

Where:
$ \phi_{xx} = \frac{\partial^2 \phi}{\partial x^2}, \quad \phi_{yy} = \frac{\partial^2 \phi}{\partial y^2}, \quad \phi_{zz} = \frac{\partial^2 \phi}{\partial z^2} $
$ \phi_{xy} = \frac{\partial^2 \phi}{\partial x \partial y}, \quad \phi_{yz} = \frac{\partial^2 \phi}{\partial y \partial z}, \quad \phi_{zx} = \frac{\partial^2 \phi}{\partial z \partial x} $

These second-order partial derivatives are also computed using central differencing.

#### Constructive Solid Geometry (CSG)

CSG is a technique used in computer graphics to create complex surfaces and solids by combining simpler ones using Boolean operations like union, intersection, and difference. Each object is represented as an implicit surface, and the operations are applied by combining their implicit functions.

#### Implementation Overview

1. **Gradient Calculation (`GetGradient`)**:

   - Uses central differencing to compute the partial derivatives of the implicit function in x, y, and z directions.
   - Constructs the gradient vector from these derivatives.
2. **Curvature Calculation (`GetCurvature`)**:

   - Computes the first-order partial derivatives using central differencing.
   - Computes the second-order partial derivatives using central differencing.
   - Uses these derivatives to calculate the mean curvature based on the formula derived from differential geometry.

# lab5

### 1.Differential Functions

A level set represents an implicit surface defined by a scalar function $\phi(x, y, z)$. The level set method is powerful for modeling and animation as it allows for robust deformations and natural topology changes. The key to level set methods is solving partial differential equations (PDEs) to evolve the surface over time.

#### 2. Differential Calculations

**Objective:** Implement various differential operators to compute gradients and other derivatives necessary for level set evolution.

1. **Differential Operators:**

   We need to implement the following differential operators:
   $$
   \phi_x^+, \phi_x^-, \phi_x^{\pm}, \phi_{xx}^{\pm}, \phi_{xy}^{\pm}
   $$
   Here, $ \phi_x^+ $ and $ \phi_x^- $ denote forward and backward differences, respectively, while $ \phi_x^{\pm} $ represents the central difference.
2. **Class Structure:**

   The functions are implemented in the `LevelSet` class with the naming convention `Diff + {X,Y,Z} + {p,m,pm}` where:

   - `{X,Y,Z}` indicates the dimension along which the derivative is taken.
   - `p` stands for plus (forward difference), `m` for minus (backward difference), and `pm` for central difference.

   ```cpp
   float LevelSet::DiffXPlus(int i, int j, int k);
   float LevelSet::DiffXMinus(int i, int j, int k);
   float LevelSet::DiffXCentral(int i, int j, int k);
   ```

#### 3. Implementing Differential Functions

1. **Forward Difference (Plus):**

   ```cpp
   float LevelSet::DiffXPlus(int i, int j, int k) {
       return (mGrid[i+1][j][k] - mGrid[i][j][k]) / mDx;
   }
   ```
2. **Backward Difference (Minus):**

   ```cpp
   float LevelSet::DiffXMinus(int i, int j, int k) {
       return (mGrid[i][j][k] - mGrid[i-1][j][k]) / mDx;
   }
   ```
3. **Central Difference:**

   ```cpp
   float LevelSet::DiffXCentral(int i, int j, int k) {
       return (mGrid[i+1][j][k] - mGrid[i-1][j][k]) / (2 * mDx);
   }
   ```
4. **Second-Order Central Difference:**

   ```cpp
   float LevelSet::DiffXXCentral(int i, int j, int k) {
       return (mGrid[i+1][j][k] - 2 * mGrid[i][j][k] + mGrid[i-1][j][k]) / (mDx * mDx);
   }
   ```

#### 4. Verifying Differential Functions

1. **Gradient Verification:**

   - Convert an implicit surface to a level set object using "Convert implicit to levelset".
   - Add a "Vector cut plane" to visualize gradients.
   - Implement and verify `LevelSet::GetGradient()` using the differential functions.
2. **Reinitialization:**

   - Reinitialize the level set to ensure it maintains the signed distance property.
   - Check the gradient and level sets after reinitialization.

#### 5. Erosion and Dilation

1. **Objective:**

   - Implement erosion and dilation operators using the upwind scheme and Godunov's method for stability.
2. **Speed Function:**

   - Erosion and dilation are achieved by setting a uniform speed function $ F $:
     $
     \text{Dilation: } F(x) = c_1, \quad c_1 > 0
     $
     $
     \text{Erosion: } F(x) = c_2, \quad c_2 < 0
     $
3. **Implementation:**

   - `ComputeTimestep()` using Equation 12 for stability.
   - `Evaluate()` to return the rate of change $ \frac{\partial \phi}{\partial t} $.

   ```cpp
   float OperatorDilateErode::Evaluate(int i, int j, int k) {
       // Compute the rate of change using upwind scheme
       return LevelSetOperator::Godunov(i, j, k);
   }
   ```
4. **Verification:**

   - Add a level set object and apply the erosion and dilation operations.
   - Check the changes using cut planes and visual inspections.

---

### 2. The Signed Distance Property

A signed distance function (SDF) is a scalar field where each point’s value represents the shortest distance to the surface, with the sign indicating whether the point is inside or outside the surface. This property is crucial for level set methods as it ensures well-defined gradients and numerical stability during surface evolution.

#### 2. Visualizing the Level Set Function

1. **Add a Level Set Object:**

   - Convert an implicit surface to a level set object using the menu option “Convert implicit to levelset”.
   - This conversion samples the implicit function and stores the values in a grid.
2. **Add a Scalar Cut Plane:**

   - Add a scalar cut plane to the level set using the corresponding menu option.
   - Select the “Iso contour” colormap, which maps scalar values to colors cycling from green to red. Each green-red cycle represents a change in the level set function by 0.1.
3. **Observe the Visualization:**

   - The scalar cut plane visualizes all level sets of the function, not just the zero level set.
   - The zero level set is rendered as a surface, and it can adjust the opacity slider to reveal the interior of the surface mesh.

#### 3. Reinitializing the Level Set Function

1. **Run the Reinitialize Operator:**

   - Execute the reinitialize operator on the level set.
   - This operation adjusts the level set function to approximate a signed distance function.
2. **Effects of Reinitialization:**

   - After reinitialization, observe how the scalar cut plane and the grid lines (spaced by 0.1) in the GUI align with the level set function.
   - Proper reinitialization ensures that the gradients are well-defined and the level sets are evenly spaced.

#### 4. Importance of Reinitialization

1. **Purpose of Reinitialization:**

   - Reinitialization maintains the level set function close to a signed distance function, which is essential for accurate numerical computations.
   - It corrects distortions that occur during surface evolution, preserving numerical stability and accuracy.
2. **Expected Results:**

   - After reinitialization, gradients should be consistent, and level sets should appear smooth and evenly spaced.
   - This ensures accurate computation of geometric properties like curvature and normals.
3. **Necessity of Reinitialization:**

   - Without reinitialization, numerical errors can cause the level set function to deviate from a signed distance function.
   - This deviation can lead to inaccurate surface evolution and instability in the numerical schemes.

### Implementation Steps

1. **Convert Implicit to Levelset:**

   - Use the “Convert implicit to levelset” menu option to convert an implicit surface to a level set object.
2. **Add a Scalar Cut Plane:**

   - Select the level set object and add a scalar cut plane using the menu option.
   - Use the “Iso contour” colormap for visualization.
3. **Run Reinitialize Operator:**

   - Execute the reinitialize operator on the level set object to adjust it to a signed distance function.
4. **Verify the Results:**

   - Check the scalar cut plane and vector cut plane to ensure that the gradients and level sets are correctly reinitialized.
   - Observe how the visualization aligns with the grid lines and verify that the gradients are well-defined and level sets are smooth.

---

### 3.Advection

#### 1. **Advection Equation and Level Set Method**

The primary goal is to advect (move) a level set surface in a vector field. The level set equation used for advection is:

$ \frac{\partial \phi}{\partial t} + \mathbf{V} \cdot \nabla \phi = 0 $

where:

- $\phi$ is the level set function.
- $\mathbf{V}$ is the velocity field.
- $\nabla \phi$ is the gradient of the level set function.

This equation describes how the level set function $\phi$ evolves over time under the influence of the velocity field $\mathbf{V}$.

#### 2. **Courant-Friedrichs-Lewy (CFL) Condition**

For the numerical solution of hyperbolic PDEs like the advection equation to be stable, the time step $\Delta t$ must satisfy the CFL condition:

$ \Delta t < \min \left( \frac{\Delta x}{|V_x|}, \frac{\Delta y}{|V_y|}, \frac{\Delta z}{|V_z|} \right) $

where $\Delta x$, $\Delta y$, and $\Delta z$ are the spatial step sizes in each direction, and $|V_x|$, $|V_y|$, and $|V_z|$ are the absolute values of the velocity components.

The CFL condition ensures that the numerical wave does not travel more than one grid cell per time step, maintaining stability in the numerical solution.

#### 3. **ComputeTimestep Function**

The `ComputeTimestep` function calculates a stable time step based on the CFL condition.

- Version 1: Calculates the Euclidean norm (magnitude) of the maximum velocity vector and divides the grid spacing by this value to determine the time step.
- Version 2: Determines the maximum component-wise velocity and uses this to compute the time step, including a safety factor.

Both approaches ensure that the chosen time step respects the CFL condition, maintaining stability during advection.

#### 4. **Upwind Scheme for Gradient Computation**

The upwind scheme is used to compute the gradient of the level set function $\nabla \phi$. This scheme ensures stability by considering the direction of the velocity field:

$$
\text{Gradient Component} = \begin{cases}
\phi_{i+1} - \phi_i & \text{if } V > 0 \\
\phi_i - \phi_{i-1} & \text{if } V < 0
\end{cases}
$$

This method ensures that the gradient is computed using values from the "upwind" direction, i.e., from the direction the wave is coming from, thus avoiding numerical instabilities that could arise from using central differences in hyperbolic equations.

#### 5. **Propagate Function and Time Integration**

The `Propagate` function advances the level set function in time using the computed stable time step. It ensures the level set function evolves correctly over the specified time period using numerical integration methods:

- **Euler Integration**: A simple method that computes the next value using the current value and the rate of change.
- **Runge-Kutta Integration**: A more accurate method that uses intermediate steps to compute the next value, reducing numerical error.

The combination of these integration methods with the CFL condition ensures that the numerical solution is both stable and accurate.

#### 6. **Evaluation Function**

The `Evaluate` function computes the rate of change $\frac{\partial \phi}{\partial t}$ at each grid point by transforming grid coordinates to world coordinates, sampling the velocity field, and computing the gradient using the upwind scheme. This rate of change is used by the integration methods to update the level set function.

### Conclusion

The code works cohesively with the tasks because it adheres to the fundamental principles of numerical stability and accuracy in solving hyperbolic PDEs. By:

- Correctly computing a stable time step using the CFL condition.
- Using the upwind scheme for gradient computation.
- Applying appropriate time integration methods.

### 4. mean curvature flow

1. **Compute the First Derivatives**:

   - `dx`, `dy`, `dz` are the central differences in the x, y, and z directions, respectively, representing the first derivatives of the level set function $\phi$.
2. **Compute the Second Derivatives**:

   - `dxx`, `dyy`, `dzz` are the second-order derivatives in the x, y, and z directions.
   - `dyz`, `dzx`, `dxy` are the mixed second-order derivatives.
3. **Calculate Gradient Magnitude Squared**:

   - Compute the sum of the squares of the first derivatives: $dx^2 + dy^2 + dz^2$.
4. **Calculate the Denominator**:

   - The denominator is $2 \cdot (dx^2 + dy^2 + dz^2)^{1.5}$.
5. **Compute Curvature Terms**:

   - Calculate each term of the curvature using the second derivatives and mixed derivatives.
6. **Sum the Curvature Terms**:

   - Sum up the curvature terms to get the mean curvature $\kappa$.
7. **Compute Gradient Magnitude**:

   - Calculate the magnitude of the gradient: $\sqrt{dx^2 + dy^2 + dz^2}$.
8. **Return Rate of Change**:

   - Finally, the rate of change is given by $\alpha \cdot \kappa \cdot |\nabla \phi|$.

# lab6

### 1.  basic functionality

> Fluid Simulation

#### 1. Introduction

In fluid simulation, we aim to model the behavior of fluids like water, smoke, and fire using computational techniques. This lab focuses on implementing a basic fluid solver for the Navier-Stokes equations using the Stable Fluids approach, which ensures numerical stability and efficiency.

### 2. The Navier-Stokes Equations

The Navier-Stokes equations for incompressible flow describe the change in the velocity field $\mathbf{V}$ of a fluid over time:

$
\frac{\partial \mathbf{V}}{\partial t} = -(\mathbf{V} \cdot \nabla)\mathbf{V} + \nu \nabla^2 \mathbf{V} - \nabla p + \mathbf{F}
$
$
\nabla \cdot \mathbf{V} = 0
$

- **Self-Advection ($(\mathbf{V} \cdot \nabla)\mathbf{V}$):** This term moves the fluid with itself, creating swirls and vortices.
- **Viscous Diffusion ($\nu \nabla^2 \mathbf{V}$):** This term represents internal friction, affecting fluid thickness (high in honey, low in water).
- **Pressure Gradient ($\nabla p$):** Ensures fluid incompressibility by adjusting pressure.
- **External Forces ($\mathbf{F}$):** Forces like gravity and wind that act on the fluid.

### 3. Solving the Navier-Stokes Equations

To solve these equations, we use operator splitting to break them into simpler sub-problems:

1. **Self-Advection:**
   $
   \frac{\partial \mathbf{V}_1}{\partial t} = -(\mathbf{V} \cdot \nabla)\mathbf{V}
   $

   - **Method:** Backward tracing to solve the advection term.
   - **Equation:**
     $
     \mathbf{V}_1(x) = \mathbf{V}(T(x, -\Delta t))
     $
   - **Interpolation:** Trilinear interpolation ensures smoothness.
2. **External Forces:**
   $
   \frac{\partial \mathbf{V}_2}{\partial t} = \mathbf{F}
   $

   - **Method:** Explicit Euler integration.
   - **Equation:**
     $
     \mathbf{V}_2 = \mathbf{V}_1 + \Delta t \cdot \mathbf{F}
     $
3. **Projection:**
   $
   \nabla \cdot \mathbf{V}_3 = 0
   $

   - **Method:** Helmholtz-Hodge decomposition to enforce incompressibility.
   - **Equation:**
     $
     \mathbf{V}_3 = \mathbf{V}_2 - \nabla q
     $
   - **Poisson Equation:** Solve for the scalar field $ q $ (pressure).
     $
     \nabla^2 q = \nabla \cdot \mathbf{V}_2
     $

### 4. Mathematical Formulation

1. **Self-Advection:**

   - **Goal:** Move the velocity field along itself.
   - **Equation:**
     $
     \frac{\partial \mathbf{V}}{\partial t} = -(\mathbf{V} \cdot \nabla)\mathbf{V}
     $
   - **Backward Tracing:**
     $
     \mathbf{V}(x, t + \Delta t) = \mathbf{V}(x - \mathbf{V}(x, t) \cdot \Delta t, t)
     $
2. **External Forces:**

   - **Goal:** Apply external forces like gravity.
   - **Equation:**
     $
     \mathbf{V} = \mathbf{V} + \Delta t \cdot \mathbf{F}
     $
   - **Explicit Euler Integration:**
     $
     \mathbf{V}(x, t + \Delta t) = \mathbf{V}(x, t) + \Delta t \cdot \mathbf{F}(x, t)
     $
3. **Projection Step:**

   - **Goal:** Ensure the velocity field is divergence-free.
   - **Equation:**
     $
     \mathbf{V} = \mathbf{V} - \nabla p
     $
   - **Poisson Equation for Pressure:**
     $
     \nabla^2 p = \nabla \cdot \mathbf{V}
     $
   - **Discrete Divergence:**
     $
     \nabla \cdot \mathbf{V} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} + \frac{\partial w}{\partial z}
     $
   - **Discrete Laplacian:**
     $
     \nabla^2 p \approx \frac{p_{i+1,j,k} + p_{i-1,j,k} + p_{i,j+1,k} + p_{i,j-1,k} + p_{i,j,k+1} + p_{i,j,k-1} - 6p_{i,j,k}}{\Delta x^2}
     $

### 5. Implementation Steps

1. **Self-Advection:**

   - Implement backward tracing to move the velocity field along itself.
   - Use trilinear interpolation for smooth transitions.
2. **External Forces:**

   - Apply external forces using explicit Euler integration.
   - Ensure forces are accurately added to the velocity field.
3. **Projection:**

   - Solve the Poisson equation to compute the pressure field.
   - Subtract the gradient of the pressure field from the velocity field to ensure incompressibility.

### 2.  semi-Lagrangian

This method is crucial for simulating inviscid (zero viscosity) fluids and ensures numerical stability even with larger time steps. The semi-Lagrangian approach is a key component in the Stable Fluids method.

### 1. Purpose of Semi-Lagrangian Integration

1. **Stability:**

   - Traditional advection methods can become unstable with large time steps. The semi-Lagrangian method is unconditionally stable, allowing for larger time steps without numerical instability.
2. **Accuracy:**

   - By tracing particles backward in time, we can accurately interpolate the velocity field, preserving the detailed flow structures like vortices and swirls.
3. **Efficiency:**

   - Semi-Lagrangian methods decouple the time step size from the spatial resolution, making the simulation more efficient.

### 2. Advection in Fluid Simulation

#### 2.1. The Advection Equation

The self-advection term in the Navier-Stokes equations is given by:
$
\frac{\partial \mathbf{V}}{\partial t} = -(\mathbf{V} \cdot \nabla) \mathbf{V}
$
This term describes the movement of the velocity field along itself.

#### 2.2. Semi-Lagrangian Method

The semi-Lagrangian method involves tracing a particle backward in time from its current position to find where it came from, then interpolating the velocity at that point.

### 3. Mathematical Formulation

#### 3.1. Backward Tracing

1. **Backward Tracing:**

   - For a given grid point $\mathbf{x}$, trace backward along the velocity field to find the departure point $\mathbf{x}_d$.
     $
     \mathbf{x}_d = \mathbf{x} - \mathbf{V}(\mathbf{x}, t) \Delta t
     $
2. **Iterative Tracing:**

   - Perform multiple small steps for better accuracy:
     $
     \mathbf{x}_{d_{i+1}} = \mathbf{x}_d - \mathbf{V}(\mathbf{x}_{d_i}, t) \frac{\Delta t}{\text{steps}}
     $
   - Here, $\text{steps}$ is the number of sub-steps in the backward tracing.

#### 3.2. Velocity Interpolation

1. **Trilinear Interpolation:**
   - Once the departure point $\mathbf{x}_d$ is found, interpolate the velocity at this point using trilinear interpolation:
     $
     \mathbf{V}_{new}(\mathbf{x}) = \mathbf{V}(\mathbf{x}_d)
     $
   - This ensures smooth transitions and avoids numerical artifacts.

#### Results

1. **Numerical Viscosity:**

   - Despite modeling inviscid fluids, the semi-Lagrangian method introduces numerical viscosity due to interpolation errors, leading to slightly diffusive results.
2. **Stability:**

   - The semi-Lagrangian method remains stable even with larger time steps, avoiding the instability issues present in purely Eulerian advection schemes.
3. **Efficiency:**

   - Larger time steps reduce the number of iterations required, improving computational efficiency.

The semi-Lagrangian method is an essential tool for fluid simulation, providing a stable and efficient way to solve the self-advection term in the Navier-Stokes equations. By tracing particles backward in time and interpolating velocities, we achieve accurate and realistic fluid behavior while maintaining numerical stability. This method is particularly useful for simulating inviscid fluids, despite the presence of some numerical viscosity due to interpolation.

# Appendix

### Helmholtz-Hodge Decomposition

#### Introduction to Helmholtz-Hodge Decomposition

The Helmholtz-Hodge decomposition theorem is a fundamental result in vector calculus. It states that any sufficiently smooth vector field $\mathbf{V}$ in a bounded domain can be decomposed into a divergence-free component, a curl-free component, and a harmonic component. For incompressible fluid simulation, we focus on the divergence-free and curl-free components.

#### Mathematical Formulation

Given a vector field $\mathbf{V}$, the Helmholtz-Hodge decomposition can be expressed as:
$
\mathbf{V} = \mathbf{V}_d + \nabla q
$
where:

- $\mathbf{V}_d$ is the divergence-free component ($\nabla \cdot \mathbf{V}_d = 0$).
- $\nabla q$ is the curl-free component ($\nabla \times \nabla q = 0$), where $q$ is a scalar potential field.

#### Why Helmholtz-Hodge Decomposition?

In fluid simulation, we need to ensure that the velocity field $\mathbf{V}$ is divergence-free to satisfy the incompressibility condition ($\nabla \cdot \mathbf{V} = 0$). The Helmholtz-Hodge decomposition allows us to modify the velocity field by subtracting the gradient of a scalar field, ensuring the result is divergence-free.

#### Applying Helmholtz-Hodge Decomposition

To apply the Helmholtz-Hodge decomposition in our fluid solver, we perform the following steps:

1. **Calculate the Divergence of the Velocity Field:**
   $
   \nabla \cdot \mathbf{V} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} + \frac{\partial w}{\partial z}
   $
2. **Solve the Poisson Equation for the Scalar Potential $ q $:**
   $
   \nabla^2 q = \nabla \cdot \mathbf{V}
   $

   - The Poisson equation can be solved using iterative methods such as Conjugate Gradient.
3. **Subtract the Gradient of $ q $ from the Velocity Field:**
   $
   \mathbf{V}_d = \mathbf{V} - \nabla q
   $

   - This step ensures that the resulting velocity field $\mathbf{V}_d$ is divergence-free.

#### Discrete Formulation

In a discrete setting, we use finite differences to approximate the operators:

1. **Discrete Divergence:**
   $
   (\nabla \cdot \mathbf{V})_{i,j,k} \approx \frac{u_{i+1,j,k} - u_{i-1,j,k}}{2\Delta x} + \frac{v_{i,j+1,k} - v_{i,j-1,k}}{2\Delta y} + \frac{w_{i,j,k+1} - w_{i,j,k-1}}{2\Delta z}
   $
2. **Discrete Laplacian:**
   $
   (\nabla^2 q)_{i,j,k} \approx \frac{q_{i+1,j,k} + q_{i-1,j,k} + q_{i,j+1,k} + q_{i,j-1,k} + q_{i,j,k+1} + q_{i,j,k-1} - 6q_{i,j,k}}{\Delta x^2}
   $
3. **Solving the Discrete Poisson Equation:**

   - Use an iterative solver like Conjugate Gradient to solve the linear system:
     $
     A q = b
     $
     where $A$ represents the discrete Laplacian operator and $b$ represents the discrete divergence of the velocity field.
4. **Update the Velocity Field:**
   $
   \mathbf{V}_d = \mathbf{V} - \nabla q
   $

   - Use finite differences to compute the gradient of $q$:
     $
     (\nabla q)_{i,j,k} \approx \left( \frac{q_{i+1,j,k} - q_{i-1,j,k}}{2\Delta x}, \frac{q_{i,j+1,k} - q_{i,j-1,k}}{2\Delta y}, \frac{q_{i,j,k+1} - q_{i,j,k-1}}{2\Delta z} \right)
     $

The Helmholtz-Hodge decomposition is a powerful tool in fluid simulation, allowing us to enforce the incompressibility condition by decomposing the velocity field into a divergence-free component and a curl-free component. By solving a Poisson equation and subtracting the gradient of the scalar potential, we ensure that the fluid’s velocity field remains divergence-free, leading to stable and accurate simulations.

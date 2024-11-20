# Ecuaciones Diferenciales Parciales mediante el método de Elemento Finito.
Este repositorio consiste de desarrollo de código para el cálculo de EDP por medio del método de elemento finito apoyándonos del cómputo científico. Este código está enfocado en la implementación de elementos finitos usando polinomios de Lagrange y cuadraturas de Gauss-Legendre. La primera parte del código consiste en la construcción de matrices de rigidez y masa, pues con ellas es que resolvemos los problemas. También en esta parte se incluye la aplicación de condiciones de Dirichlet así cómo el cálculo de la norma $L$ y $H$. 

### 1. Bibliotecas Requeridas

El código comienza con la importación de tres bibliotecas fundamentales:

- **`numpy`**: Utilizada para manejar operaciones con vectores y matrices, esencial en cálculos numéricos que involucran grandes conjuntos de datos.
- **`matplotlib.pyplot`**: Es utilizada para visualizar gráficamente los resultados, especialmente las comparaciones entre soluciones analíticas y aproximadas.
- **`math`**: Se utiliza para realizar operaciones matemáticas estándar como funciones trigonométricas o exponenciales.
- **`pandas`**: Aunque no se usa en el fragmento, esta librería normalmente se usa para manejar y analizar datos.

### 2. Raíces de Gauss-Legendre

La función **`p_roots(q)`** calcula los nodos (raíces) y los pesos de la cuadratura de Gauss-Legendre para un polinomio de grado `q`. Esto es clave para la integración numérica, ya que la cuadratura de Gauss-Legendre utiliza nodos y pesos específicos para aproximar la integral de una función en el intervalo $[-1, 1]$.

- **Entrada**: `q` es el grado del polinomio, entre 1 y 7.
- **Salida**: La función devuelve dos arrays de `numpy`:
  - `xq`: Las raíces del polinomio de Legendre.
  - `wq`: Los pesos asociados a cada nodo.

El código cubre casos de `q` entre 1 y 7, cada uno con un conjunto específico de nodos y pesos predefinidos para optimizar la precisión de la cuadratura.

### 3. Cuadratura Gaussiana

La función **`gaussian_quad(f, xq, wq)`** implementa la fórmula de la cuadratura de Gauss-Legendre para aproximar la integral de una función `f` sobre el intervalo $[-1, 1]$.

La integral se aproxima de la siguiente manera:
$$\int_{-1}^{1} f(x) \, dx \approx \sum_{i=1}^{q} w_i f(x_i)$$
- **Entrada**:
  - `f`: Función a integrar.
  - `xq`: Nodos de la cuadratura de Gauss-Legendre.
  - `wq`: Pesos asociados a los nodos.
- **Salida**: La aproximación de la integral.

### 4. Transformación y su Inversa

Estas funciones se utilizan para cambiar el intervalo de integración de $[-1, 1]$ a cualquier otro intervalo $[a, b]$.

- **`T(x, a, b)`**: Realiza la transformación lineal de un punto `x` del intervalo estándar $[-1, 1]$ al intervalo $[a, b]$.
  
  La fórmula es:
  $$T(x) = \frac{(b-a)x + (b+a)}{2}$$
- **`T_inv(x, a, b)`**: Realiza la transformación inversa de un punto `x` del intervalo $[a, b]$ al intervalo estándar $[-1, 1]$.
  
  La fórmula es:
  $$T^{-1}(x) = \frac{2}{b-a} \left( x - \frac{a+b}{2} \right)$$

### 5. Polinomios de Lagrange y sus Derivadas

Los **polinomios de Lagrange** son utilizados para interpolar datos y en este caso, para construir funciones base en un método de elementos finitos. Existen diferentes polinomios de Lagrange para diferentes grados $P$.

- **Polinomios de grado 1, 2 y 3**: Se definen funciones lambda que representan estos polinomios para valores específicos de $x$.
- **Derivadas de los polinomios**: También se definen las derivadas de los polinomios de Lagrange para grados 1, 2 y 3, lo cual es importante cuando se trabaja con la matriz de rigidez, ya que implica el uso de derivadas de los polinomios de base.

### 6. Diccionarios de Polinomios y Derivadas

Se crean dos diccionarios:
- **`MD`**: Almacena los polinomios de Lagrange para diferentes grados.
- **`dMD`**: Almacena las derivadas de los polinomios correspondientes.

Estos diccionarios permiten una consulta eficiente de los polinomios y sus derivadas según el grado y el índice.

### 7. Evaluación de Polinomios de Lagrange y Derivadas

- **`L(x, l, p)`**: Esta función evalúa el polinomio de Lagrange de grado $p$ en el punto $x$ para el índice $l$. Consulta el diccionario `MD` para obtener el polinomio correspondiente y luego lo evalúa.
  
- **`dL(x, l, p)`**: Similar a la función anterior, pero evalúa la derivada del polinomio de Lagrange. Consulta el diccionario `dMD` para obtener la derivada y luego la evalúa.

### 8. Matrices de Masa y Rigidez Locales

Las matrices de masa y rigidez son componentes fundamentales en los métodos de elementos finitos. Estas matrices se construyen a partir de las funciones base (polinomios de Lagrange) y sus derivadas.

- **Matriz de Masa (`M_loc`)**: Representa la integral de los productos de los polinomios de Lagrange:
  $$M_{ij} = \int_{-1}^{1} L_i(x) \cdot L_j(x) \, dx$$
  Esta matriz se arma usando la cuadratura de Gauss-Legendre para calcular las integrales.
  
- **Matriz de Rigidez (`S_loc`)**: Similar a la de masa, pero involucra las derivadas de los polinomios de Lagrange:
  $$S_{ij} = \int_{-1}^{1} dL_i(x) \cdot dL_j(x) \, dx$$
  
Estas matrices se utilizan en los métodos de elementos finitos para representar la relación entre las fuerzas y los desplazamientos en los problemas de mecánica estructural o física computacional.

Claro, a continuación te doy una descripción detallada de las funciones que mencionas:

### 9. **RHS (Right Hand Side)**

El vector del lado derecho (RHS) en problemas de elementos finitos generalmente representa la parte del problema que no depende de la solución directa, como las fuentes o términos externos que se deben considerar. En este contexto, el RHS generalmente incluye la aplicación de una función externa en el dominio y su integración sobre los elementos. La función `RHS` calcula el vector del lado derecho para un sistema de ecuaciones discretizado. Dependiendo del problema, puede estar relacionada con un término de fuente $f(x)$, es decir, la integral de $f(x)\cdot L_i(x)$, donde $L_i(x)$ es un polinomio de Lagrange.

**Entrada:**
- `f : callable`: Función a integrar, la cual puede representar una fuente o término en la ecuación.
- `p : int`: Grado del polinomio de Lagrange.
- `xq : numpy.ndarray`: Nodos de cuadratura de Gauss.
- `wq : numpy.ndarray`: Pesos de la cuadratura de Gauss.

**Salida:**
- `rhs : numpy.ndarray`: Vector del lado derecho que contiene las integrales para cada uno de los elementos.

### 10. **Aplicar Condiciones de Dirichlet**

Las condiciones de Dirichlet se aplican cuando conocemos el valor de la solución en los bordes del dominio. Esto se realiza estableciendo las correspondientes entradas en el sistema de ecuaciones a cero y asignando los valores conocidos a las incógnitas. La función `apply_dirichlet_conditions` aplica las condiciones de Dirichlet a un sistema de ecuaciones. Esto implica modificar la matriz global de rigidez y el vector del lado derecho para que los nodos de frontera tengan valores prescritos.

**Entrada:**
- `K : numpy.ndarray`: Matriz global de rigidez.
- `rhs : numpy.ndarray`: Vector global del lado derecho.
- `dirichlet_nodes : list`: Índices de los nodos donde se aplican las condiciones de Dirichlet.
- `dirichlet_values : list`: Valores prescritos en esos nodos.

**Salida:**
- `K, rhs`: Matriz de rigidez modificada y el vector del lado derecho con las condiciones de Dirichlet aplicadas.

### 11. **Norma $L^2$**

La norma $L^2$ se utiliza para medir el "tamaño" de una función en el espacio $L^2$, que es la norma asociada a la integral del cuadrado de la función.

#### Descripción
La función `L2_Norm` calcula la norma $ L^2 $ de una función aproximada utilizando los polinomios de Lagrange. La norma se define como:

$$\| u \|_{L^2} = \left( \int_a^b |u(x)|^2 dx \right)^{1/2}$$

**Entrada:**
- `u : callable`: La función aproximada cuya norma $ L^2 $ se quiere calcular.
- `xq : numpy.ndarray`: Nodos de cuadratura de Gauss.
- `wq : numpy.ndarray`: Pesos de la cuadratura de Gauss.

**Salida:**
- `norm_L2 : float`: El valor de la norma $ L^2 $.


### 12. **Norma $ H^1 $**

La norma $H^1$ es una norma del espacio de Sobolev, que es una extensión de la norma $L^2$ que también toma en cuenta las derivadas de la función.

#### Descripción
La función `H_Norm` calcula la norma $ H^1 $, que se define como:

$$\| u \|_{H^1} = \left( \| u \|_{L^2}^2 + \| u' \|_{L^2}^2 \right)^{1/2}$$

Esto mide no solo la magnitud de la función, sino también la de su derivada.

**Entrada:**
- `u : callable`: La función aproximada cuya norma $H^1$ se quiere calcular.
- `du : callable`: La derivada de la función $u$.
- `xq : numpy.ndarray`: Nodos de cuadratura de Gauss.
- `wq : numpy.ndarray`: Pesos de la cuadratura de Gauss.

**Salida:**
- `norm_H : float`: El valor de la norma $H^1$.


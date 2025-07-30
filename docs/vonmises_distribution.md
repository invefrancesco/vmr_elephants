

# Distribuzione Von Mises

## Introduzione

Il post di [Whittenbury
Daniel](https://dlwhittenbury.github.io/ds-1-sampling-and-visualising-the-von-mises-distribution.html)
contiene alcune note riguardanti la distribuzione von Mises.

La “statistica direzionale” ha numerosi usi nelle scienze

- La statistica classica, come sviluppata da Gauss, assume che gli
  errori di misura siano “piccoli” e distribuiti in uno **spazio lineare
  infinito** (come la retta dei numeri reali o uno spazio euclideo)
- Questo funziona bene quando le misure sono **molto precise** come nel
  caso di astronomia o dei geometria che usano strumenti accurati
- **Le misure angolari** (come direzioni su una bussola, ore del giorno,
  orientamenti spaziali, ecc.) **vivono naturalmente su una struttura
  topologica diversa** (un cerchio o una sfera).
- Se gli errori diventano grandi (cioè, se i dati sono molto dispersi),
  **non ha più senso ignorare la curvatura** dello spazio su cui si
  trovano i dati.
- Per esempio, dire che la media tra 10° e 350° è 180° è assurdo in un
  contesto circolare, ma potrebbe sembrare corretto se si ignorasse la
  topologia.

La distribuzione di von Mises è stata introdotta da Richard von Mises
nel 1918 e la sua PDF è

- x è la variabile aleatoria angolare
- $\mu$ è la **direzione media** (con $\mu \in [0, 2 \pi)$)
- $\kappa$ è un parametro di concentrazione della distribuzione, che può
  assumere valori non negativi ($\kappa \in [0, \infty)$)
  - più è grande $\kappa$, più la distribuzione è concentrata attorno
    alla direzione media
- $I_0(\kappa)$ è la funzione di Bessel modificata del primo tipo e
  ordine zero

Le distribuzioni circolari hanno diversi vantaggi

- sono definite su uno spazio campionario corretto
  - **la normale può produrre stime errate vicino ai limiti (es. 0 e
    $2\pi$)**
- la distribuzione di von Mises ha un comportamento limite corretto
  quando $\kappa \to 0$: tende a una distribuzione uniforme sul cerchio
  dove tutte le direzioni sono equiprobabili
- **usare la normale in problemi circolari può portare a errori seri
  come fondere in modo scorretto una distribuzione a priori con una
  funzione di verosimiglianza**
  - **nelle misurazioni dell’angolo di arrivo (in applicazioni di
    tracciamento), si rischia di aggiornare male la stima bayesiana se
    si ignora la natura ciclica della variabile**

## Optim function

A common task in statistics is maximizing (or minimizing) complex
univariate or multivariate functions. This is typically done in the
context of maximizing likelihood functions with respect to a vector of
unknown parameters, given the observed data.

- More generally, this task pertains to an “optimization” problem, where
  we seek to finding a set of parameter values that minimizes or
  maximizes some pre-defined objective function.
- When this objective function is a log likelihood function, this
  reduces to the problem of maximum likelihood estimation.
- There are many potential approaches to use for optimization in this
  context.
- How best to optimize a function of interest depends on the nature of
  the function to be optimized, as well as practical concerns regarding
  a candidate procedure to be used.

By default `optim` performs minimization, but it will maximize if
control\$fnscale is negative. optimHess is an auxiliary function to
compute the Hessian at a later stage if hessian = TRUE was forgotten.

- The default method is an implementation of that of Nelder and Mead
  (1965), that uses only function values and is robust but relatively
  slow. It will work reasonably well for non-differentiable functions.
- Method “BFGS” is a quasi-Newton method (also known as a variable
  metric algorithm), specifically that published simultaneously in 1970
  by Broyden, Fletcher, Goldfarb and Shanno. This uses function values
  and gradients to build up a picture of the surface to be optimized.
- Method “CG” is a conjugate gradients method based on that by Fletcher
  and Reeves (1964) (but with the option of Polak–Ribiere or
  Beale–Sorenson updates). Conjugate gradient methods will generally be
  more fragile than the BFGS method, but as they do not store a matrix
  they may be successful in much larger optimization problems.
- Method “L-BFGS-B” is that of Byrd et al. (1995) which allows box
  constraints, that is each variable can be given a lower and/or upper
  bound. The initial value must satisfy the constraints. This uses a
  limited-memory modification of the BFGS quasi-Newton method. If
  non-trivial bounds are supplied, this method will be selected, with a
  warning.
- Nocedal and Wright (1999) is a comprehensive reference for the
  previous three methods.
- Method “SANN” is by default a variant of simulated annealing given in
  Belisle (1992).
- Simulated-annealing belongs to the class of stochastic global
  optimization methods. It uses only function values but is relatively
  slow.
- It will also work for non-differentiable functions. This
  implementation uses the Metropolis function for the acceptance
  probability.
- By default the next candidate point is generated from a Gaussian
  Markov kernel with scale proportional to the actual temperature. If a
  function to generate a new candidate point is given, method “SANN” can
  also be used to solve combinatorial optimization problems. -
  Temperatures are decreased according to the logarithmic cooling
  schedule as given in Belisle (1992, p. 890); specifically, the
  temperature is set to `temp / log(((t-1) %/% tmax) * tmax + exp(1))`,
  where `t` is the current iteration step and `temp` and `tmax` are
  specifiable via control, see below. - Note that the “SANN” method
  depends critically on the settings of the control parameters. It is
  not a general-purpose method but can be very useful in getting to a
  good value on a very rough surface. - Method “Brent” is for
  one-dimensional problems only, using
  `optimize(<ff>, lower, upper, tol = control$reltol)` where `<ff>` is
  `function(par) fn(par, ...)/control$fnscale`. It can be useful in
  cases when `optim()` is used inside other functions where only method
  can be specified, such as in `mle` from package stats4.

Function `fn` can return `NA` or `Inf` if the function cannot be
evaluated at the supplied value, but the initial value must have a
computable finite value of fn. (Except for method “L-BFGS-B” where the
values should always be finite.)

`optim` can be used recursively, and for a single parameter as well as
many. It also accepts a zero-length par, and just evaluates the function
with that argument

## Stima di una Von Mises con una covariata

Creo due vettori $\mathbf{x}_n$ e $\mathbf{y}_n$ in cui

- $X_i \sim N(0, 1)$
- $Y_i \sim VonMises(\mu = 2atan(\beta_0 +\beta_1 \cdot x_i), k = 10)$
  con $\beta_0=-1$, $\beta_1 = 3.5$

Definisco la funzione di - log-verosimiglianza

Minimizzo la funione con `optim`

    [1] TRUE

| Parameter | Estimate | Std. Error | 95% CI Lower | 95% CI Upper |
|:----------|:--------:|:----------:|:------------:|:------------:|
| beta 0    |  -1.002  |   0.019    |    -1.040    |    -0.964    |
| beta 1    |  3.519   |   0.054    |    3.413     |    3.624     |
| kappa     |  2.315   |   0.043    |    2.229     |    2.400     |

Nota esplicativa:

- Il parametro `kappa` è stimato in scala logaritmica per garantire che
  il valore restituito sia sempre positivo. Il valore `par[3]`
  rappresenta quindi il logaritmo naturale del parametro reale, e viene
  esponenziato all’interno della funzione di verosimiglianza.
- Impostando `hessian = TRUE`, `optim()` restituisce anche la **matrice
  Hessiana** della funzione obiettivo valutata nel punto di minimo, che
  corrisponde (sotto opportune condizioni regolarità) alla **matrice di
  informazione osservata**.
- Se questa matrice Hessiana è invertibile, la sua inversa fornisce una
  **stima asintotica della matrice di varianza-covarianza** dei
  parametri stimati.
- Le **radici quadrate degli elementi diagonali** della matrice inversa
  corrispondono agli **errori standard asintotici** delle stime dei
  parametri.

## Simulazione di $K$ repliche per lo stesso set di parametri veri

Genero $K$ set di dati random e stimo i parametri che massimizzano la
verosimiglianza per ogni set di dati

Per riproducibilità genero una sola volta una lista di $K = 1000$
`tibble` inserendo come come parametri `n = 1000`, `beta.0 = -1`,
`beta.1 = 3.5`, `kappa = 10` e la salvo in un file chiamato
`sim_1.RData`

| Parameter |  Mean   | Std. Dev |  2.5%   |  97.5%  |
|:----------|:-------:|:--------:|:-------:|:-------:|
| beta 0    | -1.0001 |  0.0198  | -1.0386 | -0.9616 |
| beta 1    | 3.5010  |  0.0558  | 3.3933  | 3.6093  |
| kappa     | 10.0305 |  0.4347  | 9.2221  | 10.9423 |

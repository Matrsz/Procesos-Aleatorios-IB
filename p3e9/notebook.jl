### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 451411e9-6eb1-4ad1-97cd-946300511c0d
begin
	import Pkg
	Pkg.add(url="https://github.com/Pocket-titan/DarkMode")
    import DarkMode
end

# ╔═╡ 53d9d2c0-0369-403e-8378-858663608a35
begin
    Pkg.activate(mktempdir())
    Pkg.add([
		Pkg.PackageSpec(name="Revise"),
        Pkg.PackageSpec(name="Latexify"),
        Pkg.PackageSpec(name="LaTeXStrings"),
        Pkg.PackageSpec(name="Plots"),
		Pkg.PackageSpec(name="ColorSchemes"),
		Pkg.PackageSpec(name="PlotThemes"),
		Pkg.PackageSpec(name="Statistics")
    ])
	using Revise
    using Latexify
    using LaTeXStrings
	using Plots
	using ColorSchemes
	using PlotThemes
	using Statistics
end

# ╔═╡ d4e2266c-9103-43ce-b54c-1885f12c8297
md"## P3E9 - Procesos Estocásticos

Considere una señal de amplitud modulada $v(t) = \sum_{n=-\infty}^\infty I_n g(t − nT)$ donde $I_n$ corresponde a un proceso ESA de valores discretos con función de autocorrelación $R_I[m]$ y $g(t)$ es una función determinística cualquiera.

En particular tomamos:

$$\begin{align}
	\bullet\quad &I_n \sim \mathcal U {[-1,\,1]}\hspace{27em}\\[1em]
	\bullet\quad &g(t)=\Pi(t-\tfrac 1 2)\\[1em]
	\bullet\quad &T = 1
\end{align}$$"

# ╔═╡ c09b4e07-1d9b-45d6-8a1b-ca6a1321d0e0
begin
	step(t) = t >= 0 ? 1 : 0
	box(t) = step(t+1/2) - step(t-1/2)
	g(t) = box(t-1/2)
end

# ╔═╡ 8840e441-7553-413a-aa00-7a431533fda0
I(n)=rand(-1:1)

# ╔═╡ a4cc31e0-77af-41d0-8cb4-4ada252a119b
T = 1

# ╔═╡ 32dc34a7-d85e-4aeb-8bc0-fb14f5c79e06
md"#### Grafique una realización del proceso.
Definimos en función de $I_n,\; g(t),\; T$ genéricos una función que es capaz de retornar una instancia genérica de $v(t)$, y la usamos para instanciar el caso particular de $v(t)$ con los $I_n,\; g(t),\; T$ dados"

# ╔═╡ 4e9d3f16-1728-4e59-9ff0-b8116da9a7e8
function specify_v(I, g, T)
	n = -100:100
	function v_instance(t) 
		A = zeros(size(t))
		for i in n
			A .+= I(i)*g.(t.-i*T)
		end
		return A
	end
	return v_instance
end

# ╔═╡ 1e7309eb-fc00-4898-9047-e90ee61a6255
v = specify_v(I, g, T)

# ╔═╡ 9b938890-ac7a-467d-889b-b7dd4644f175
md"Usando esta instancia de $v(t)$, simulamos algunas realizaciones." 

# ╔═╡ 5c775a41-0b9e-4dc9-85a6-b50cf8461c33
md"#### Obtenga analíticamente la autocorrelación del proceso $v(t)$

Partiendo de la definición genérica de $v(t)$

$\begin{aligned}
R_v(t_1,t_2) &\textstyle= E\left[v(t_1)\cdot v(t_2)\right]\\[1em]
R_v(t_1,t_2) &\textstyle = E\left[\sum_{n=-\infty}^\infty I_n g(t_1 − nT)\cdot\sum_{m=-\infty}^\infty I_m g(t_2 − mT)\right]\\[1em]
R_v(t_1,t_2) &\textstyle = E\left[\sum_{n=-\infty}^\infty\sum_{m=-\infty}^\infty I_n I_m  g(t_1 − nT) g(t_2 − mT)\right]\\[1em]
R_v(t_1,t_2) &\textstyle = \sum_{n=-\infty}^\infty\sum_{m=-\infty}^\infty E\left[I_n I_m  g(t_1 − nT) g(t_2 − mT)\right]\\[1em]
R_v(t_1,t_2) &\textstyle = \sum_{n=-\infty}^\infty\sum_{m=-\infty}^\infty E\left[I_n I_m\right]  g(t_1 − nT) g(t_2 − mT)\\[1em]
R_v(t_1,t_2) &\textstyle = \sum_{n=-\infty}^\infty\sum_{m=-\infty}^\infty R_I[n,m]  g(t_1 − nT) g(t_2 − mT)
\end{aligned}$

A continuación se especializa $I_n$, considerando que es

+ Independiente: $E[I_nI_m] = E[I_n]E[I_m]$
+ Idénticamente distribuido: $E[I_n] = \mu, \; Var[I_n] = \sigma^2$
$\begin{aligned}
R_I[n,m] &= \left\lvert\;
			\begin{aligned}
				\rlap{E[I_nI_m]}\hphantom{Var[I_n]+E[I_n]^2} \qquad &si \quad n\ne m\\
				\rlap{E[I_n^2]}\hphantom{Var[I_n]+E[I_n]^2} \qquad &si \quad n=m
\end{aligned}\right.\\[0.5em]
R_I[n,m] &= \left\lvert\;
			\begin{aligned}
				\rlap{E[I_n]E[I_m]}\hphantom{Var[I_n]+E[I_n]^2} \qquad &si \quad n\ne m\\
				Var[I_n]+E[I_n]^2 \qquad &si \quad n=m
\end{aligned}\right.\\[0.5em]
R_I[n,m] &= \left\lvert\;
			\begin{aligned}
				\rlap{\mu^2}\hphantom{Var[I_n]+E[I_n]^2} \qquad &si \quad n\ne m\\
				\rlap{\sigma^2+\mu^2}\hphantom{Var[I_n]+E[I_n]^2} \qquad &si \quad n=m
\end{aligned}\right.\\[0.5em]
R_I[n,m] &= \mu^2 + \delta_{nm}\sigma^2
\end{aligned}$
Ya que $I_n\sim\mathcal U{[-1,1]}$ se conocen $\mu = \frac{b-a}2=0, \; \sigma^2 = \frac{n^2-1}{12}=\frac 2 3$. Resulta $R_I[n,m] = \frac 2 3 \delta_{mn}$

$
\begin{aligned}
R_v(t_1,t_2)  &\textstyle = \sum_{n=-\infty}^\infty\sum_{m=-\infty}^\infty \frac 2 3 \delta_{nm}  g(t_1 − nT) g(t_2 − mT)\\[1em]
R_v(t_1,t_2)  &\textstyle = \frac 2 3 \sum_{n=-\infty}^\infty g(t_1 − nT) g(t_2 − nT)
\end{aligned}$

Se especializa finalmente $g(t) = \Pi(t-\tfrac 1 2)$ usando la definición de la función rectangular en 2 dimensiones, $\Pi(x,y)=\Pi(x)\Pi(y)$, para llegar al resultado final:

$\boxed{R_v(t_1,t_2) = \frac 2 3 \sum_{n=-\infty}^\infty \Pi\left(t_1- nT-\tfrac 1 2,t_2- nT-\tfrac 1 2\right)}$
"

# ╔═╡ bb8ef4f5-b8f6-41e6-b623-0c50559c3f0f
md"#### Determine si el proceso es cicloestacionario

Un proceso $X(t)$ se dice cicloestacionario en sentido amplio con periodo $T_0$ si se cumplen

$\begin{align}
	\bullet\quad &E[X(t+T_0)] = E[X(t)]\hspace{27em}\\[1em]
	\bullet\quad &R_X[t_1+T_0, t_2+T_0] = R_X(t_1,t_2)
\end{align}$

Para $v(t)$, el primer inciso se cumple trivialmente, ya que vimos anteriormente que $E[v(t)] = 0$

Para evaluar la segunda propiedad asumiendo cicloestacionariedad de periodo $T$ usamos

$\begin{aligned}
R_v(t_1+T,t_2+T) &= \frac 2 3 \sum_{n=-\infty}^\infty g(t_1 + T − nT) g(t_2 + T − nT)\\[0.5em]
R_v(t_1+T,t_2+T) &= \frac 2 3 \sum_{n=-\infty}^\infty g(t_1 − (n-1)T) g(t_2 − (n-1)T)\\[0.5em]
R_v(t_1+T,t_2+T) &= \frac 2 3 \sum_{k=-\infty}^\infty g(t_1 − kT) g(t_2 − kT) \qquad con \quad k=n-1\\[1em]
R_v(t_1+T,t_2+T) &= R_v(t_1,t_2) 
\end{aligned}$

Por ende, el proceso $v(t)$ es cicloestacionario en sentido amplio, independiente de la defnición de $g(t)$
"

# ╔═╡ 50c77bf9-4c69-4dc5-bf22-f057f2d502eb
md"#### Obtenga la función de autocorrelación en forma empírica

Para calcular la autocorrelación empírica definimos la función autocorr, que calcula la autocorrelación empírica de una función aleatoria genérica *f* para tiempos *t* con *nReal* realizaciones:"

# ╔═╡ bf850a6e-9a1f-4aec-9f1d-f832df3779dc
function autocorr(f, t, nReal)
	N = length(t)
	R = zeros(N, N)     
	for i in 1:nReal
		y = f(t)       
		R += y*y'/nReal		 
	end
	return R
end

# ╔═╡ cba35460-eb6e-4c29-9052-4dd5630055fd
md"El resultado es consistente con el resultado analítico obtenido en el inciso anterior"

# ╔═╡ 14f2c153-4b24-4697-b6a8-6e9b1ac2bd33
md"#### Se desconoce el instante inicial, verificar que el proceso es ESA

Si desconocemos el instante inicial en que la señal fué generada podemos modelar a la señal agregando un retardo aleatorio, obteniendo $v_2(t) = \sum_{n=-\infty}^\infty I_n g(t − nT+\phi)$ con $\phi \sim \mathcal U{(0, T)}$

Podemos definir un nuevo proceso aleatorio $f(t) = g(t+\phi) \longrightarrow v_2(t) = \sum_{n=-\infty}^\infty I_n f(t − nT)$

$\begin{aligned}
R_{v_2}(t_1,t_2) &\textstyle= E\left[v_2(t_1)\cdot v_2(t_2)\right]\\[1em]
R_{v_2}(t_1,t_2) &\textstyle = E\left[\sum_{n=-\infty}^\infty I_n f(t_1 − nT)\cdot\sum_{m=-\infty}^\infty I_m f(t_2 − mT)\right]\\[1em]
R_{v_2}(t_1,t_2) &\textstyle = \sum_{n=-\infty}^\infty\sum_{m=-\infty}^\infty E\left[I_n I_m  f(t_1 − nT) f(t_2 − mT)\right]\\[1em]
R_{v_2}(t_1,t_2) &\textstyle = \sum_{n=-\infty}^\infty\sum_{m=-\infty}^\infty E\left[I_n I_m\right]E\left[  f(t_1 − nT) f(t_2 − mT)\right]\\[1em]
R_{v_2}(t_1,t_2) &\textstyle = \sum_{n=-\infty}^\infty\sum_{m=-\infty}^\infty R_I[m.n]R_f(t_1-nT, t_2-mT)\\
&\vdots \\
R_{v_2}(t_1,t_2)  &\textstyle = \frac 2 3 \sum_{n=-\infty}^\infty R_f(t_1-nT, t_2-nT)
\end{aligned}$

Encontramos la autocorrelación de $f(t)$

$\begin{aligned}
R_{f}(t_1,t_2) &\textstyle= E\left[g(t_1+\phi)\cdot g(t_2+\phi)\right]\\[1em]
R_{f}(t_1,t_2) &\textstyle = \int_\infty^\infty g(t_1+\phi) g(t_2+\phi)\,p(\phi)\,d\phi\\[1em]
R_{f}(t_1,t_2) &\textstyle = \frac 1 T \int_0^T g(t_1+\phi) g(t_2+\phi)\,d\phi
\end{aligned}$

La autocorrelación de $v_2(t)$ resulta entonces

$\textstyle R_{v_2}(t_1,t_2)  = \frac 2 3 \sum_{n=-\infty}^\infty \frac 1 T \int_0^T g(t_1-nT+\phi) g(t_2-nT+\phi)\,d\phi$

Aplicando cambio de variables: $\phi = \xi - t_2,\; \tau = t_1-t_2$

$\begin{aligned}
R_{v_2}(t_1,t_2)  &\textstyle = \frac 2 3 \frac 1 T \sum_n \int_{t_2}^{T+t_2} g(\tau-nT+\xi) g(\xi-nT)\,d\xi\\[1em]
R_{v_2}(t_1,t_2)  &\textstyle = \small \frac 2 3 \frac 1 T \left[\sum_n \int_{t_2}^{T} g(\tau-nT+\xi) g(\xi-nT)\,d\xi + \int_{T}^{T+t_2} g(\tau-nT+\xi) g(\xi-nT)\,d\xi \right]\\[1em] 
R_{v_2}(t_1,t_2)  &\textstyle = \small \frac 2 3 \frac 1 T \left[\sum_{n}\int_{t_2}^{T} g(\tau-nT+\xi) g(\xi-nT)\,d\xi + \sum_{n} \int_{T}^{T+t_2} g(\tau-nT+\xi) g(\xi-nT)\,d\xi \right]\\[1em]
R_{v_2}(t_1,t_2)  &\textstyle = \small \frac 2 3 \frac 1 T \left[\sum_{n} \int_{t_2}^{T} g(\tau-nT+\xi) g(\xi-nT)\,d\xi + \sum_{n} \int_{0}^{t_2} g(\tau-nT+\gamma-T) g(\gamma-T-nT)\,d\gamma\right] \\[1em]
R_{v_2}(t_1,t_2)  &\textstyle = \small \frac 2 3 \frac 1 T \left[\sum_{n} \int_{t_2}^{T} g(\tau-nT+\xi) g(\xi-nT)\,d\xi +  \sum_{k} \int_{0}^{t_2} g(\tau-kT+\gamma) g(\gamma-kT)\,d\gamma\right]\\[1em]
R_{v_2}(t_1,t_2)  &\textstyle = \small \frac 2 3 \frac 1 T \left[\sum_{n} \int_{t_2}^{T} g(\tau-nT+\xi) g(\xi-nT)\,d\xi + \int_{0}^{t_2} g(\tau-kT+\gamma) g(\gamma-kT)\,d\gamma\right]\\[1em]
R_{v_2}(t_1,t_2)  &\textstyle = \frac 2 3 \frac 1 T \sum_{n} \int_{0}^{T} g(\tau-nT+\xi) g(\xi-nT)\,d\xi\\[1em]
R_{v_2}(t_1,t_2)  &\textstyle = R_{v_2}(\tau, 0)
\end{aligned}$

Demostrando que $v_2(t)$ resulta siempre estacionario en sentido amplio, independientemente de la definición de $g(t)$
"

# ╔═╡ 9760d2c8-59e4-4dba-9d69-d98c13d02dac
md"#### Obtenga la autocorrelación de $v_2(t)$ en forma empírica

Para simular $v_2(t)$ solo hace falta modificar ligeramente la función utilizada para instanciar $v(t)$ definida anteriormente, agregando el término $\phi$"

# ╔═╡ ef3f534f-102e-482d-a805-673719658095
function specify_v2(I, g, T)
	n = -100:100
	function v_instance(t) 
		ϕ = rand(0:0.001:T)
		A = zeros(size(t))
		for i in n
			A .+= I(i)*g.(t.-i*T.+ϕ)
		end
		return A
	end
	return v_instance
end

# ╔═╡ 2f90fced-c538-4eac-85b1-238fc3b85c89
v2 = specify_v2(I, g, T)

# ╔═╡ 8b12d927-a056-4987-a654-0bfd15af6a99
md"Instanciamos algunas realizaciones de la función $v_2(t)$"

# ╔═╡ d66656d5-34d5-43d4-b03f-20ebc21c1aa8
md"Y finalmente encontramos la autocorrelación empírica de $v_2(t)$, usando la misma función *autocorr* definida anteriormente"

# ╔═╡ b69f9cbc-b116-4206-9ede-2a3a2d6bc135
begin
	struct Wow
		filename
	end

	function Base.show(io::IO, ::MIME"image/png", w::Wow)
		write(io, read(w.filename))
	end
end

# ╔═╡ f7f0cca8-e2d0-42a4-baf3-06cf132585f5
begin
	plt = :seaborn_deep6
	cgr = :acton
	DarkMode.enable(theme="shadowfox")
end

# ╔═╡ ffbce74c-6791-416e-a474-4abaebbd4803
begin
	x = 0:0.001:50
	theme(:juno)
	p1 = plot(x, v(x), palette=plt, yticks=-1:1)
	p2 = plot(x, v(x), palette=plt, yticks=-1:1)
	p3 = plot(x, v(x), palette=plt, yticks=-1:1)
	p4 = plot(x, v(x), palette=plt, yticks=-1:1)
	plot(p1, p2, p3, p4, layout=(4,1), legend=false, size=(800,500))
end

# ╔═╡ c717d020-2894-41d5-b032-3741f6f5c4b1
begin
	t = 0:0.1:15
	Rv = autocorr(v, t, 5000)
	theme(:juno)
	heatmap(t, t, Rv, c=cgr, size=(800,600), title="Rv", xlabel="t1", ylabel="t2")
end

# ╔═╡ 49f0cf7f-52f3-4f38-8185-0cb898fa6b4b
begin
	F=25
	x2 = -2:0.001:F
	theme(:juno)
	p1_2 = plot(x2, v2(x2), palette=plt, yticks=-1:1, xlims=(0,20))
	p2_2 = plot(x2, v2(x2), palette=plt, yticks=-1:1, xlims=(0,20))
	p3_2 = plot(x2, v2(x2), palette=plt, yticks=-1:1, xlims=(0,20))
	p4_2 = plot(x2, v2(x2), palette=plt, yticks=-1:1, xlims=(0,20))
	plot(p1_2, p2_2, p3_2, p4_2, layout=(4,1), legend=false, size=(800,500))
end

# ╔═╡ e8a89ff1-74d5-4107-b3e8-7bac5197a83e
begin
	Rv2 = autocorr(v2, t, 5000)
	theme(:juno)
	heatmap(t, t, Rv2, c=cgr, size=(800,600), title="Rv2", xlabel="t1", ylabel="t2")
end

# ╔═╡ Cell order:
# ╟─d4e2266c-9103-43ce-b54c-1885f12c8297
# ╠═c09b4e07-1d9b-45d6-8a1b-ca6a1321d0e0
# ╠═8840e441-7553-413a-aa00-7a431533fda0
# ╟─a4cc31e0-77af-41d0-8cb4-4ada252a119b
# ╟─32dc34a7-d85e-4aeb-8bc0-fb14f5c79e06
# ╠═4e9d3f16-1728-4e59-9ff0-b8116da9a7e8
# ╠═1e7309eb-fc00-4898-9047-e90ee61a6255
# ╟─9b938890-ac7a-467d-889b-b7dd4644f175
# ╟─ffbce74c-6791-416e-a474-4abaebbd4803
# ╟─5c775a41-0b9e-4dc9-85a6-b50cf8461c33
# ╟─bb8ef4f5-b8f6-41e6-b623-0c50559c3f0f
# ╟─50c77bf9-4c69-4dc5-bf22-f057f2d502eb
# ╠═bf850a6e-9a1f-4aec-9f1d-f832df3779dc
# ╟─c717d020-2894-41d5-b032-3741f6f5c4b1
# ╟─cba35460-eb6e-4c29-9052-4dd5630055fd
# ╟─14f2c153-4b24-4697-b6a8-6e9b1ac2bd33
# ╟─9760d2c8-59e4-4dba-9d69-d98c13d02dac
# ╠═ef3f534f-102e-482d-a805-673719658095
# ╠═2f90fced-c538-4eac-85b1-238fc3b85c89
# ╟─8b12d927-a056-4987-a654-0bfd15af6a99
# ╟─49f0cf7f-52f3-4f38-8185-0cb898fa6b4b
# ╟─d66656d5-34d5-43d4-b03f-20ebc21c1aa8
# ╟─e8a89ff1-74d5-4107-b3e8-7bac5197a83e
# ╟─b69f9cbc-b116-4206-9ede-2a3a2d6bc135
# ╟─53d9d2c0-0369-403e-8378-858663608a35
# ╟─451411e9-6eb1-4ad1-97cd-946300511c0d
# ╟─f7f0cca8-e2d0-42a4-baf3-06cf132585f5

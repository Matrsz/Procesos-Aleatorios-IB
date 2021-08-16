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
        Pkg.PackageSpec(name="ControlSystems"),
		Pkg.PackageSpec(name="ColorSchemes"),
		Pkg.PackageSpec(name="PlotThemes"),
		Pkg.PackageSpec(name="Distributions"),
		Pkg.PackageSpec(name="DSP"),
		Pkg.PackageSpec(name="Statistics")
    ])
	using Revise
    using Latexify
    using LaTeXStrings
	using Plots
	using ControlSystems
	using ColorSchemes
	using PlotThemes
	using Distributions
	using DSP
	using Statistics
end

# ╔═╡ b9d08d19-ab2c-4ed8-b9c3-4cfd774a7875
L"Y_n=-\sum_{i=1}^q \alpha_i Y_{n-i} + \sum_{j=0}^p \beta_j W_{n-j}"

# ╔═╡ d6fcf049-c911-4194-9bf6-51e5afd46de1
md"Donde $W \sim N(1,0)$, que para la simulación se realiza de la siguiente manera:"

# ╔═╡ e629d4ad-e056-4895-abc5-f1b4dc153d85
begin
	d = Normal(0,1)
	W = rand(d, 50)
end

# ╔═╡ 3ab1a1d9-376b-4596-939b-886197e2ad7e
md"Para implementar la simulación en un lenguaje de computación científica, resulta util vectorizar la expresión del proceso ARMA."

# ╔═╡ c9a885da-3f5f-469a-8302-e207d7a3da08
L"Y_n=
\begin{bmatrix}
\alpha_{1} & \alpha_{2} & \cdots & \alpha_{q}
\end{bmatrix}
\begin{bmatrix}
Y_{n-1} \\ Y_{n-2} \\ \cdots \\ Y_{n-q}
\end{bmatrix}
+
\begin{bmatrix}
\beta_{0} & \beta_{1} & \cdots & \beta_{p}
\end{bmatrix}
\begin{bmatrix}
W_{n} \\ W_{n-1} \\ \cdots \\ W_{n-p}
\end{bmatrix}
"

# ╔═╡ fde9f15f-9253-4e36-9bb9-db0108860f77
function ARMA(α, β, W)
    q = length(α)
	p = length(β)
	N = length(W)
	Y = copy(W)
    for n in 1:N
        Y[n] = (n-q>0 ? α'*Y[n-1:-1:n-q] : 0) + (n-p+1>0 ? β'*W[n:-1:n-p+1] : 0)
    end
    return Y
end

# ╔═╡ ec60d196-0361-4d4d-abb5-e5ff0f095e6e
md"### Proceso AR(1)
Tomamos el caso particular proceso autoregresivo de orden 1: $q=1,\, p=0$"

# ╔═╡ 288d9139-6f59-425f-b6d7-83296239c48a
L"Y_n=
\alpha_1 Y_{n-1}
+
\beta_0 W_0
"

# ╔═╡ 057f086f-d255-4c12-b17a-b3c7849405ba
begin
	Y1 = ARMA([0.1], [1], W)
	Y2 = ARMA([0.5], [1], W)
	Y3 = ARMA([0.9], [1], W)
	Y4 = ARMA([-0.5], [1], W)
end

# ╔═╡ 9bf15602-5b55-4728-bec4-8ed2435ef274
md"Para encontrar la autocorrelación empírica de los procesos, la siguiente función retorna la matriz autocorrelación del proceso, realizando **nReal** realizaciones del mismo y promediando su valor."

# ╔═╡ 0c6b0755-c441-4bc5-b403-e45b4a908ed3
function ARMA_corr(α, β, N, nReal)
	R = zeros(N, N, nReal)      #Instancio matrix 3d de N x N x nReal
	for i in 1:nReal
		W = rand(d, N)			#Vector W de ruido blanco
		Y = ARMA(α, β, W)       #Proceso ARMA usando ese W
		R[:, :, i] = Y*Y'		#Guardo la realización de Y*Y'
	end
	Ryy = mean!(ones(N, N), R)  #Promedio las matrices (sintaxis depende el lenguaje)
	return Ryy
end

# ╔═╡ d473ef0a-9474-49ae-96be-df786f88f9cb
begin
	N = 80
	nReal = 1000
	Ryy1 = ARMA_corr([0.1], [1], N, nReal)
	Ryy2 = ARMA_corr([0.5], [1], N, nReal)
	Ryy3 = ARMA_corr([0.9], [1], N, nReal)
	Ryy4 = ARMA_corr([-0.5], [1], N, nReal)
end

# ╔═╡ 0a5f9d33-ce43-45b5-8b8f-7430868c5f76
md"Para estudiar la transferencia del proceso consideramos la definición del proceso ARMA en el dominio de la transformada Z"

# ╔═╡ 82c0d8ca-1a8d-45b0-af14-2c5e700c180e
L"Y(z)=-\sum_{i=1}^q \alpha_i z^{-i}Y(z) + \sum_{j=0}^p \beta_j z^{-j}W(z)"

# ╔═╡ dffe59fc-858c-4681-80f3-ee0c4a3e462a
L"Y(z)\left[1+\sum_{i=1}^q \alpha_i z^{-i}\right] = W(z) \sum_{j=0}^p \beta_j z^{-j}"

# ╔═╡ a21d68a7-104b-454d-9b15-1d524a07235e
md"De esa forma se obtiene la transferencia dada por $H(z) = \frac{Y(z)}{W(z)}$"

# ╔═╡ 0acf8aea-5640-43b6-ae77-5815a53ad57d
L"H(z) = \frac{\sum_{j=0}^p \beta_j z^{-j}}{1+\sum_{i=1}^q \alpha_i z^{-i}}"

# ╔═╡ bd84ce02-faa6-4c20-94b8-6be124f103ad
md"En el caso particular de $\alpha=\begin{bmatrix}1 & -0.5 & 0.25\end{bmatrix}$ y $\beta=1$, por ejemplo, la transferencia resultaría"

# ╔═╡ 95c76b56-69f3-4d46-964f-3509fc875a23
latexify(:(H(z)=1/(1+z^-1-0.5z^-2+0.25z^-3)))

# ╔═╡ 72fa3402-1501-4c8b-b90a-f185ccae0cca
begin
	α5 = [1, -0.5, 0.25]
	β5 = [1]
	H5 = tf([β5; zeros(length(α5)+1-length(β5))], [1; α5],1)
end

# ╔═╡ 0dd331de-759c-436a-b0a6-0445abef0225
pzmap(H)

# ╔═╡ df1a9787-4cc3-4a70-b36e-586d00190f69
begin
	w = range(-pi, stop=pi, length=200)
	Hw = freqz(PolynomialRatio(β, [1; α]), w)
	plot(w, abs.(Hw))
end

# ╔═╡ 4fe9cc74-a00c-4ea8-9d0d-db11d17734fb
md"### Procesos de media movil
Repetimos el proceso considerando procesos de media móvil, $q = 0,\, p \ne 0$"

# ╔═╡ e3ab60ef-0357-451b-8b2f-f388123a54fe
begin
	β9  = [1, 1]
	β10 = [1, -1]
	β11 = [0.8, 0.4, 0.2, 0.1]
	β12 = [0.8, -0.4, 0.2, -0.1]
	Y9  = ARMA([0], β9 , W)
	Y10 = ARMA([0], β10, W)
	Y11 = ARMA([0], β11, W)
	Y12 = ARMA([0], β12, W)
end

# ╔═╡ 96467507-b69a-4698-a7ae-cbb99d441a38
begin
	Ryy9  = ARMA_corr([0], β9 , N, nReal)
	Ryy10 = ARMA_corr([0], β10, N, nReal)
	Ryy11 = ARMA_corr([0], β11, N, nReal)
	Ryy12 = ARMA_corr([0], β12, N, nReal)
end

# ╔═╡ b69f9cbc-b116-4206-9ede-2a3a2d6bc135
begin
	struct Wow
		filename
	end

	function Base.show(io::IO, ::MIME"image/png", w::Wow)
		write(io, read(w.filename))
	end
end

# ╔═╡ f99891ef-ab72-4bcb-8e14-10a70cdc5750
DarkMode.enable(theme="shadowfox")

# ╔═╡ 4530598e-5fe7-4e7e-99ca-ff1adb1a6040
begin
	plt = :seaborn_deep6
	cgr = :acton
end

# ╔═╡ b3111ac1-1838-460d-a94d-f4dd6567aa2c
begin 
	theme(:juno)
	p1 = plot([W Y1], label=["W" "Y"], title="α=0.1", palette=plt)
	p2 = plot([W Y2], label=["W" "Y"], title="α=0.5", palette=plt)
	p3 = plot([W Y3], label=["W" "Y"], title="α=0.9", palette=plt)
	p4 = plot([W Y4], label=["W" "Y"], title="α=-0.5", palette=plt)
	plot(p1, p2, p3, p4, layout=(2,2))
end

# ╔═╡ 9ae4c725-4681-4b8f-8572-8bd1854e882d
begin
	theme(:juno)
	pR1 = heatmap(1:N, 1:N, Ryy1, colorbar=false, title="α=0.1", c=cgr)
	pR2 = heatmap(1:N, 1:N, Ryy2, colorbar=false, title="α=0.5", c=cgr)
	pR3 = heatmap(1:N, 1:N, Ryy3, colorbar=false, title="α=0.9", c=cgr)
	pR4 = heatmap(1:N, 1:N, Ryy4, colorbar=false, title="α=-0.5", c=cgr)
	plot(pR1, pR2, pR3, pR4, layout=(2,2))
end

# ╔═╡ fb0d90ce-9b20-45ce-98f7-f79014e0224c
begin 
	theme(:juno)
	p9  = plot([W Y9 ], label=["W" "Y"], title="β="*string(β9 ), palette=plt)
	p10 = plot([W Y10], label=["W" "Y"], title="β="*string(β10), palette=plt)
	p11 = plot([W Y11], label=["W" "Y"], title="β="*string(β11), palette=plt)
	p12 = plot([W Y12], label=["W" "Y"], title="β="*string(β12), palette=plt)
	plot(p9, p10, p11, p12, layout=(2,2))
end

# ╔═╡ 12fb86d5-a41f-445b-8254-555f34629c12
begin
	theme(:juno)
	pR9  = heatmap(1:N, 1:N, Ryy9 , colorbar=false, title="β="*string(β9 ), c=cgr)
	pR10 = heatmap(1:N, 1:N, Ryy10, colorbar=false, title="β="*string(β10), c=cgr)
	pR11 = heatmap(1:N, 1:N, Ryy11, colorbar=false, title="β="*string(β11), c=cgr)
	pR12 = heatmap(1:N, 1:N, Ryy12, colorbar=false, title="β="*string(β12), c=cgr)
	plot(pR9, pR10, pR11, pR12, layout=(2,2))
end

# ╔═╡ Cell order:
# ╟─b9d08d19-ab2c-4ed8-b9c3-4cfd774a7875
# ╟─d6fcf049-c911-4194-9bf6-51e5afd46de1
# ╠═e629d4ad-e056-4895-abc5-f1b4dc153d85
# ╟─3ab1a1d9-376b-4596-939b-886197e2ad7e
# ╟─c9a885da-3f5f-469a-8302-e207d7a3da08
# ╠═fde9f15f-9253-4e36-9bb9-db0108860f77
# ╠═ec60d196-0361-4d4d-abb5-e5ff0f095e6e
# ╟─288d9139-6f59-425f-b6d7-83296239c48a
# ╠═057f086f-d255-4c12-b17a-b3c7849405ba
# ╠═b3111ac1-1838-460d-a94d-f4dd6567aa2c
# ╟─9bf15602-5b55-4728-bec4-8ed2435ef274
# ╠═0c6b0755-c441-4bc5-b403-e45b4a908ed3
# ╠═d473ef0a-9474-49ae-96be-df786f88f9cb
# ╠═9ae4c725-4681-4b8f-8572-8bd1854e882d
# ╟─0a5f9d33-ce43-45b5-8b8f-7430868c5f76
# ╟─82c0d8ca-1a8d-45b0-af14-2c5e700c180e
# ╟─dffe59fc-858c-4681-80f3-ee0c4a3e462a
# ╟─a21d68a7-104b-454d-9b15-1d524a07235e
# ╠═0acf8aea-5640-43b6-ae77-5815a53ad57d
# ╠═bd84ce02-faa6-4c20-94b8-6be124f103ad
# ╠═95c76b56-69f3-4d46-964f-3509fc875a23
# ╠═72fa3402-1501-4c8b-b90a-f185ccae0cca
# ╠═0dd331de-759c-436a-b0a6-0445abef0225
# ╠═df1a9787-4cc3-4a70-b36e-586d00190f69
# ╠═4fe9cc74-a00c-4ea8-9d0d-db11d17734fb
# ╠═e3ab60ef-0357-451b-8b2f-f388123a54fe
# ╠═fb0d90ce-9b20-45ce-98f7-f79014e0224c
# ╠═96467507-b69a-4698-a7ae-cbb99d441a38
# ╠═12fb86d5-a41f-445b-8254-555f34629c12
# ╠═b69f9cbc-b116-4206-9ede-2a3a2d6bc135
# ╠═451411e9-6eb1-4ad1-97cd-946300511c0d
# ╠═f99891ef-ab72-4bcb-8e14-10a70cdc5750
# ╠═53d9d2c0-0369-403e-8378-858663608a35
# ╠═4530598e-5fe7-4e7e-99ca-ff1adb1a6040

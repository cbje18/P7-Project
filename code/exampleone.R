spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), mean.model = list(armaOrder = c(0,1), include.mean = FALSE))
fit <- ugarchfit(spec, data, out.sample = 100, solver = "hybrid")
forecast <- ugarchforecast(fit, data, n.ahead = 10)
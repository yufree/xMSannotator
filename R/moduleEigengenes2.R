moduleEigengenes2 <-
function (expr, colors, impute = TRUE, nPC = 1, align = "along average",
    excludeGrey = FALSE, grey = ifelse(is.numeric(colors), 0, 
        "grey"), subHubs = TRUE, trapErrors = FALSE, returnValidOnly = trapErrors, 
    softPower = 6, scale = TRUE, verbose = 0, indent = 0) 
{
    spaces = indentSpaces(indent)

    if (verbose == 1) 
        printFlush(paste(spaces, "moduleEigengenes: Calculating", 
            nlevels(as.factor(colors)), "module eigengenes in given set."))
    if (is.null(expr)) {
        stop("moduleEigengenes: Error: expr is NULL. ")
    }
    if (is.null(colors)) {
        stop("moduleEigengenes: Error: colors is NULL. ")
    }
    if (is.null(dim(expr)) || length(dim(expr)) != 2) 
        stop("moduleEigengenes: Error: expr must be two-dimensional.")
    if (dim(expr)[2] != length(colors)) 
        stop("moduleEigengenes: Error: ncol(expr) and length(colors) must be equal (one color per gene).")
    if (is.factor(colors)) {
        nl = nlevels(colors)
        nlDrop = nlevels(colors[, drop = TRUE])
        if (nl > nlDrop) 
            stop(paste("Argument 'colors' contains unused levels (empty modules). ", 
                "Use colors[, drop=TRUE] to get rid of them."))
    }
    if (softPower < 0) 
        stop("softPower must be non-negative")
    alignRecognizedValues = c("", "along average")
    if (!is.element(align, alignRecognizedValues)) {
        printFlush(paste("ModulePrincipalComponents: Error:", 
            "parameter align has an unrecognised value:", align, 
            "; Recognized values are ", alignRecognizedValues))
        stop()
    }
    maxVarExplained = 10
    if (nPC > maxVarExplained) 
        warning(paste("Given nPC is too large. Will use value", 
            maxVarExplained))
    nVarExplained = min(nPC, maxVarExplained)
    modlevels = levels(factor(colors))
    if (excludeGrey) 
        if (sum(as.character(modlevels) != as.character(grey)) > 
            0) {
            modlevels = modlevels[as.character(modlevels) != 
                as.character(grey)]
        }
        else {
            stop(paste("Color levels are empty. Possible reason: the only color is grey", 
                "and grey module is excluded from the calculation."))
        }
    PrinComps = data.frame(matrix(NA, nrow = dim(expr)[[1]], 
        ncol = length(modlevels)))
    averExpr = data.frame(matrix(NA, nrow = dim(expr)[[1]], ncol = length(modlevels)))
    varExpl = data.frame(matrix(NA, nrow = nVarExplained, ncol = length(modlevels)))
    validMEs = rep(TRUE, length(modlevels))
    validAEs = rep(FALSE, length(modlevels))
    isPC = rep(TRUE, length(modlevels))
    isHub = rep(FALSE, length(modlevels))
    validColors = colors
    names(PrinComps) = paste(moduleColor.getMEprefix(), modlevels, 
        sep = "")
    names(averExpr) = paste("AE", modlevels, sep = "")
    for (i in c(1:length(modlevels))) {
        if (verbose > 1) 
            printFlush(paste(spaces, "moduleEigengenes : Working on ME for module", 
                modlevels[i]))
        modulename = modlevels[i]
        restrict1 = as.character(colors) == as.character(modulename)
        if (verbose > 2) 
            printFlush(paste(spaces, " ...", sum(restrict1), 
                "genes"))
        datModule = as.matrix(t(expr[, restrict1]))
        n = dim(datModule)[1]
        p = dim(datModule)[2]
	return(datModule)
        pc = try({
            if (nrow(datModule) > 1 && impute) {
                seedSaved = FALSE
                if (exists(".Random.seed")) {
                  saved.seed = .Random.seed
                  seedSaved = TRUE
                }
                if (any(is.na(datModule))) {
                  if (verbose > 5) 
                    printFlush(paste(spaces, " ...imputing missing data"))
                  datModule = impute.knn(datModule, k = min(10, 
                    nrow(datModule) - 1))
                  try({
                    if (!is.null(datModule$data)) 
                      datModule = datModule$data
                  }, silent = TRUE)
                }
                if (seedSaved) 
                  .Random.seed <<- saved.seed
            }
            if (verbose > 5) 
                printFlush(paste(spaces, " ...scaling"))
            if (scale) 
                datModule = t(scale(t(datModule)))
            if (verbose > 5) 
                printFlush(paste(spaces, " ...calculating SVD"))
            svd1 = svd(datModule, nu = min(n, p, nPC), nv = min(n, 
                p, nPC))
            if (verbose > 5) 
                printFlush(paste(spaces, " ...calculating PVE"))
            veMat = cor(svd1$v[, c(1:min(n, p, nVarExplained))], 
                t(datModule), use = "p")
            varExpl[c(1:min(n, p, nVarExplained)), i] = apply(veMat^2, 
                1, mean, na.rm = TRUE)
            svd1$v[, 1]
        }, silent = TRUE)
        if (class(pc) == "try-error") {
            if ((!subHubs) && (!trapErrors)) 
                stop(pc)
            if (subHubs) {
                if (verbose > 0) {
                  printFlush(paste(spaces, " ..principal component calculation for module", 
                    modulename, "failed with the following error:"))
                  printFlush(paste(spaces, "     ", pc, spaces, 
                    " ..hub genes will be used instead of principal components."))
                }
                isPC[i] = FALSE
                pc = try({
                  scaledExpr = scale(t(datModule))
                  covEx = cov(scaledExpr, use = "p")
                  modAdj = abs(covEx)^softPower
                  kIM = (apply(modAdj, 1, sum, na.rm = TRUE))^3
                  if (max(kIM, na.rm = TRUE) > 1) 
                    kIM = kIM - 1
                  kIM[is.na(kIM)] = 0
                  hub = which.max(kIM)
                  alignSign = sign(covEx[, hub])
                  alignSign[is.na(alignSign)] = 0
                  isHub[i] = TRUE
                  pcxMat = scaledExpr * matrix(kIM * alignSign, 
                    nrow = nrow(scaledExpr), ncol = ncol(scaledExpr), 
                    byrow = TRUE)/sum(kIM)
                  pcx = apply(pcxMat, 1, sum, na.rm = TRUE)
                  varExpl[1, i] = mean(cor(pcx, t(datModule), 
                    use = "p")^2, na.rm = TRUE)
                  pcx
                }, silent = TRUE)
            }
        }
        if (class(pc) == "try-error") {
            if (!trapErrors) 
                stop(pc)
            if (verbose > 0) {
                printFlush(paste(spaces, " ..ME calculation of module", 
                  modulename, "failed with the following error:"))
                printFlush(paste(spaces, "     ", pc, spaces, 
                  " ..the offending module has been removed."))
            }
            warning(paste("Eigengene calculation of module", 
                modulename, "failed with the following error \n     ", 
                pc, "The offending module has been removed.\n"))
            validMEs[i] = FALSE
            isPC[i] = FALSE
            isHub[i] = FALSE
            validColors[restrict1] = grey
        }
        else {
            PrinComps[, i] = pc
            ae = try({
		 return(datModule)
                if (isPC[i]) 
		   
                  scaledExpr = scale(t(datModule))
                averExpr[, i] = apply(scaledExpr, 1, mean, na.rm = TRUE)
                if (align == "along average") {
                  if (verbose > 4) 
                    printFlush(paste(spaces, " .. aligning module eigengene with average expression."))
		
		#return(averExpr)
                  if (cor(averExpr[, i], PrinComps[, i], use = "p") < 
                    0) 
                    PrinComps[, i] = -PrinComps[, i]
                }
                0
            }, silent = TRUE)
            if (class(ae) == "try-error") {
                if (!trapErrors) 
                  stop(ae)
                if (verbose > 0) {
                  printFlush(paste(spaces, " ..Average expression calculation of module", 
                    modulename, "failed with the following error:"))
                  printFlush(paste(spaces, "     ", ae, spaces, 
                    " ..the returned average expression vector will be invalid."))
                }
                warning(paste("Average expression calculation of module", 
                  modulename, "failed with the following error \n     ", 
                  ae, "The returned average expression vector will be invalid.\n"))
            }
            validAEs[i] = !(class(ae) == "try-error")
        }
    }
    allOK = (sum(!validMEs) == 0)
    if (returnValidOnly && sum(!validMEs) > 0) {
        PrinComps = PrinComps[, validMEs]
        averExpr = averExpr[, validMEs]
        varExpl = varExpl[, validMEs]
        validMEs = rep(TRUE, times = ncol(PrinComps))
        isPC = isPC[validMEs]
        isHub = isHub[validMEs]
        validAEs = validAEs[validMEs]
    }
    allPC = (sum(!isPC) == 0)
    allAEOK = (sum(!validAEs) == 0)
    list(eigengenes = PrinComps, averageExpr = averExpr, varExplained = varExpl, 
        nPC = nPC, validMEs = validMEs, validColors = validColors, 
        allOK = allOK, allPC = allPC, isPC = isPC, isHub = isHub, 
        validAEs = validAEs, allAEOK = allAEOK)
}

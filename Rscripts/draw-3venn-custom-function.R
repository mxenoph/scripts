function (area1, area2, area3, n12, n23, n13, n123, category = rep("",
    3), rotation = 1, reverse = FALSE, euler.d = TRUE, scaled = TRUE,
    lwd = rep(2, 3), lty = rep("solid", 3), col = rep("black",
        3), fill = NULL, alpha = rep(0.5, 3), label.col = rep("black",
        7), cex = rep(1, 7), fontface = rep("plain", 7), fontfamily = rep("serif",
        7), cat.pos = c(-40, 40, 180), cat.dist = c(0.05, 0.05,
        0.025), cat.col = rep("black", 3), cat.cex = rep(1, 3),
    cat.fontface = rep("plain", 3), cat.fontfamily = rep("serif",
        3), cat.just = list(c(0.5, 1), c(0.5, 1), c(0.5, 0)),
    cat.default.pos = "outer", cat.prompts = FALSE, rotation.degree = 0,
    rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist = 0.05,
    offset = 0, inter.lab= NULL, ...)
{
    if (length(category) == 1) {
        cat <- rep(category, 3)
    }
    else if (length(category) != 3) {
        stop("Unexpected parameter length for 'category'")
    }
    if (length(lwd) == 1) {
        lwd <- rep(lwd, 3)
    }
    else if (length(lwd) != 3) {
        stop("Unexpected parameter length for 'lwd'")
    }
    if (length(lty) == 1) {
        lty <- rep(lty, 3)
    }
    else if (length(lty) != 3) {
        stop("Unexpected parameter length for 'lty'")
    }
    if (length(col) == 1) {
        col <- rep(col, 3)
    }
    else if (length(col) != 3) {
        stop("Unexpected parameter length for 'col'")
    }
    if (length(label.col) == 1) {
        label.col <- rep(label.col, 7)
    }
    else if (length(label.col) != 7) {
        stop("Unexpected parameter length for 'label.col'")
    }
    if (length(cex) == 1) {
        cex <- rep(cex, 7)
    }
    else if (length(cex) != 7) {
        stop("Unexpected parameter length for 'cex'")
    }
    if (length(fontface) == 1) {
        fontface <- rep(fontface, 7)
    }
    else if (length(fontface) != 7) {
        stop("Unexpected parameter length for 'fontface'")
    }
    if (length(fontfamily) == 1) {
        fontfamily <- rep(fontfamily, 7)
    }
    else if (length(fontfamily) != 7) {
        stop("Unexpected parameter length for 'fontfamily'")
    }
    if (length(fill) == 1) {
        fill <- rep(fill, 3)
    }
    else if (length(fill) != 3 & length(fill) != 0) {
        stop("Unexpected parameter length for 'fill'")
    }
    if (length(alpha) == 1) {
        alpha <- rep(alpha, 3)
    }
    else if (length(alpha) != 3 & length(alpha) != 0) {
        stop("Unexpected parameter length for 'alpha'")
    }
    if (length(cat.pos) == 1) {
        cat.pos <- rep(cat.pos, 3)
    }
    else if (length(cat.pos) != 3) {
        stop("Unexpected parameter length for 'cat.pos'")
    }
    if (length(cat.dist) == 1) {
        cat.dist <- rep(cat.dist, 3)
    }
    else if (length(cat.dist) != 3) {
        stop("Unexpected parameter length for 'cat.dist'")
    }
    if (length(cat.col) == 1) {
        cat.col <- rep(cat.col, 3)
    }
    else if (length(cat.col) != 3) {
        stop("Unexpected parameter length for 'cat.col'")
    }
    if (length(cat.cex) == 1) {
        cat.cex <- rep(cat.cex, 3)
    }
    else if (length(cat.cex) != 3) {
        stop("Unexpected parameter length for 'cat.cex'")
    }
    if (length(cat.fontface) == 1) {
        cat.fontface <- rep(cat.fontface, 3)
    }
    else if (length(cat.fontface) != 3) {
        stop("Unexpected parameter length for 'cat.fontface'")
    }
    if (length(cat.fontfamily) == 1) {
        cat.fontfamily <- rep(cat.fontfamily, 3)
    }
    else if (length(cat.fontfamily) != 3) {
        stop("Unexpected parameter length for 'cat.fontfamily'")
    }
    if (!(class(cat.just) == "list" & length(cat.just) == 3)) {
        stop("Unexpected parameter format for 'cat.just'")
    }
    else if (!(length(cat.just[[1]]) == 2 & length(cat.just[[2]]) ==
        2 & length(cat.just[[3]]) == 2)) {
        stop("Unexpected parameter format for 'cat.just'")
    }
    if (euler.d == FALSE & scaled == TRUE) {
        stop("Uninterpretable parameter combination\nPlease set both euler.d = FALSE and scaled = FALSE to force Venn diagrams.")
    }
    if (offset > 1 | offset < 0) {
        stop("'Offset' must be between 0 and 1.  Try using 'rotation.degree = 180' to achieve offsets in the opposite direction.")
    }
    cat.pos <- cat.pos + rotation.degree
    a1 <- area1 - n12 - n13 + n123
    a2 <- n12 - n123
    a3 <- area2 - n12 - n23 + n123
    a4 <- n13 - n123
    a5 <- n123
    a6 <- n23 - n123
    a7 <- area3 - n13 - n23 + n123
    areas <- c(a1, a2, a3, a4, a5, a6, a7)
    if (euler.d) {
        special.code <- VennDiagram::decide.special.case(areas)
        function.name <- paste("draw.", special.code, sep = "")
        if (function.name %in% ls("package:VennDiagram")) {
            f1 <- get(function.name)
            rst <- f1(a1 = areas[1], a2 = areas[2], a3 = areas[3],
                a4 = areas[4], a5 = areas[5], a6 = areas[6],
                a7 = areas[7], category = category, reverse = reverse,
                cat.default.pos = cat.default.pos, lwd = lwd,
                lty = lty, col = col, label.col = label.col,
                cex = cex, fontface = fontface, fontfamily = fontfamily,
                cat.pos = cat.pos, cat.dist = cat.dist, cat.col = cat.col,
                cat.cex = cat.cex, cat.fontface = cat.fontface,
                cat.fontfamily = cat.fontfamily, cat.just = cat.just,
                cat.prompts = cat.prompts, fill = fill, alpha = alpha,
                ...)
            rst <- VennDiagram::adjust.venn(VennDiagram::rotate.venn.degrees(gList1 = rst,
                angle = rotation.degree, x.centre = rotation.centre[1],
                y.centre = rotation.centre[2]), ...)
            if (ind) {
                grid.draw(rst)
            }
            return(rst)
        }
    }
    rotated <- VennDiagram::rotate(areas, category, lwd, lty,
        col, label.col, cex, fontface, fontfamily, cat.col, cat.cex,
        cat.fontface, cat.fontfamily, alpha, rotation, reverse,
        fill)
    for (i in 1:length(areas)) {
        areas[i] <- rotated[[1]][i]
    }
    category <- rotated[[2]]
    lwd <- rotated$lwd
    lty <- rotated$lty
    col <- rotated$col
    label.col <- rotated$label.col
    cex <- rotated$cex
    fontface <- rotated$fontface
    fontfamily <- rotated$fontfamily
    cat.col <- rotated$cat.col
    cat.cex <- rotated$cat.cex
    cat.fontface <- rotated$cat.fontface
    cat.fontfamily <- rotated$cat.fontfamily
    fill <- rotated$fill
    alpha <- rotated$alpha
    areas.error <- c("a1 <- area1 - n12 - n13 + n123", "a2 <- n12 - n123",
        "a3 <- area2 - n12 - n23 + n123", "a4 <- n13 - n123",
        "a5 <- n123", "a6 <- n23 - n123", "a7 <- area3 - n13 - n23 + n123")
    for (i in 1:length(areas)) {
        if (areas[i] < 0) {
            stop(paste("Impossible:", areas.error[i], "produces negative area"))
        }
    }
    for (i in 1:length(areas)) {
        if (areas[i]) {
            scaled <- FALSE
        }
    }
    is.defaults <- TRUE
    if (is.expression(category)) {
        is.defaults <- FALSE
    }
    if (all(cat.default.pos != "outer", cat.default.pos != "text",
        !is.defaults, cat.prompts)) {
        print("No default location recognized.  Automatically changing to 'outer'")
        cat.default.pos <- "outer"
    }
    if (all(cat.default.pos == "outer", !is.defaults, cat.prompts)) {
        print("Placing category labels at default outer locations.  Use 'cat.pos' and 'cat.dist' to modify location.")
        print(paste("Current 'cat.pos':", cat.pos[1], "degrees,",
            cat.pos[2], "degrees"))
        print(paste("Current 'cat.dist':", cat.dist[1], ",",
            cat.dist[2]))
    }
    if (all(cat.default.pos == "text", !is.defaults, cat.prompts)) {
        print("Placing category labels at default text locations.  Use 'cat.pos' and 'cat.dist' to modify location.")
        print(paste("Current 'cat.pos':", cat.pos[1], "degrees,",
            cat.pos[2], "degrees"))
        print(paste("Current 'cat.dist':", cat.dist[1], ",",
            cat.dist[2]))
    }
    grob.list <- gList()
    if (!exists("overrideTriple")) {
        r1 <- sqrt(100/pi)
        r2 <- r1
        r3 <- r1
    }
    else {
        r1 <- sqrt(area1/pi)
        r2 <- sqrt(area2/pi)
        r3 <- sqrt(area3/pi)
    }
    max.circle.size = 0.2
    shrink.factor <- max.circle.size/r1
    r1 <- r1 * shrink.factor
    r2 <- r2 * shrink.factor
    r3 <- r3 * shrink.factor
    if (!exists("overrideTriple")) {
        a <- find.dist(100, 100, 40) * shrink.factor
        b <- a
        c <- a
    }
    else {
        a <- find.dist(area1, area2, n12) * shrink.factor
        b <- find.dist(area2, area3, n23) * shrink.factor
        c <- find.dist(area1, area3, n13) * shrink.factor
    }
    x.centres <- vector(mode = "numeric", length = 3)
    y.centres <- vector(mode = "numeric", length = 3)
    beta <- (a^2 + c^2 - b^2)/(2 * a * c)
    gamma <- sqrt(1 - beta^2)
    x.centres[1] <- (r1 - r2 - a + 1)/2
    x.centres[3] <- x.centres[1] + c * beta
    y.centres[3] <- (r3 - r1 + 1 - c * gamma)/2
    y.centres[1] <- y.centres[3] + c * gamma
    x.centres[2] <- x.centres[1] + a
    y.centres[2] <- y.centres[1]
    radii <- c(r1, r2, r3)
    for (i in 1:3) {
        grob.list <- gList(grob.list, VennDiagram::ellipse(x = x.centres[i],
            y = y.centres[i], a = radii[i], b = radii[i], gp = gpar(lty = 0,
                fill = fill[i], alpha = alpha[i])))
    }
    for (i in 1:3) {
        grob.list <- gList(grob.list, VennDiagram::ellipse(x = x.centres[i],
            y = y.centres[i], a = radii[i], b = radii[i], gp = gpar(lwd = lwd[i],
                lty = lty[i], col = col[i], fill = "transparent")))
    }
    new.x.centres <- vector(mode = "numeric", length = 3)
    new.y.centres <- vector(mode = "numeric", length = 3)
    cell.labels <- paste0(areas, inter.lab)
    cell.x <- vector(mode = "numeric", length = 7)
    cell.y <- vector(mode = "numeric", length = 7)
    x.cept.12 <- (r1^2 - r2^2 - x.centres[1]^2 + x.centres[2]^2)/(2 *
        (x.centres[2] - x.centres[1]))
    y.cept.12.1 <- sqrt(r1^2 - (x.cept.12 - x.centres[1])^2) +
        y.centres[1]
    y.cept.12.2 <- -sqrt(r1^2 - (x.cept.12 - x.centres[1])^2) +
        y.centres[1]
    theta <- acos((a^2 + c^2 - b^2)/(2 * a * c))
    new.x.centres[3] <- x.centres[1] + c
    l.x.cept.13 <- (r1^2 - r3^2 - x.centres[1]^2 + new.x.centres[3]^2)/(2 *
        (new.x.centres[3] - x.centres[1]))
    l.y.cept.13.1 <- sqrt(r1^2 - (l.x.cept.13 - x.centres[1])^2) +
        y.centres[1]
    l.y.cept.13.2 <- -sqrt(r1^2 - (l.x.cept.13 - x.centres[1])^2) +
        y.centres[1]
    rot <- sqrt(2 * r1^2 - 2 * r1^2 * cos(theta))
    x.cept.13.1 <- l.x.cept.13 + rot * cos(pi/2 - atan((l.y.cept.13.1 -
        y.centres[1])/(l.x.cept.13 - x.centres[1])) + theta/2)
    x.cept.13.2 <- l.x.cept.13 + rot * cos(pi/2 - atan((l.y.cept.13.2 -
        y.centres[1])/(l.x.cept.13 - x.centres[1])) + theta/2)
    y.cept.13.1 <- l.y.cept.13.1 - rot * sin(pi/2 - atan((l.y.cept.13.1 -
        y.centres[1])/(l.x.cept.13 - x.centres[1])) + theta/2)
    y.cept.13.2 <- l.y.cept.13.2 - rot * sin(pi/2 - atan((l.y.cept.13.2 -
        y.centres[1])/(l.x.cept.13 - x.centres[1])) + theta/2)
    theta <- -acos((a^2 + b^2 - c^2)/(2 * a * b))
    new.x.centres[3] <- x.centres[2] - b
    l.x.cept.23 <- (r2^2 - r3^2 - x.centres[2]^2 + new.x.centres[3]^2)/(2 *
        (new.x.centres[3] - x.centres[2]))
    l.y.cept.23.1 <- sqrt(r2^2 - (l.x.cept.23 - x.centres[2])^2) +
        y.centres[2]
    l.y.cept.23.2 <- -sqrt(r2^2 - (l.x.cept.23 - x.centres[2])^2) +
        y.centres[2]
    rot <- sqrt(2 * r2^2 - 2 * r2^2 * cos(theta))
    x.cept.23.1 <- l.x.cept.23 + rot * cos(pi/2 - atan((y.centres[2] -
        l.y.cept.23.1)/(x.centres[2] - l.x.cept.23)) + theta/2)
    x.cept.23.2 <- l.x.cept.23 + rot * cos(pi/2 - atan((y.centres[2] -
        l.y.cept.23.2)/(x.centres[2] - l.x.cept.23)) + theta/2)
    y.cept.23.1 <- l.y.cept.23.1 - rot * sin(pi/2 - atan((y.centres[2] -
        l.y.cept.23.1)/(x.centres[2] - l.x.cept.23)) + theta/2)
    y.cept.23.2 <- l.y.cept.23.2 - rot * sin(pi/2 - atan((y.centres[2] -
        l.y.cept.23.2)/(x.centres[2] - l.x.cept.23)) + theta/2)
    m <- (y.cept.23.2 - y.cept.23.1)/(x.cept.23.2 - x.cept.23.1)
    y.sect <- m * (x.cept.12 - x.cept.23.1) + y.cept.23.1
    cell.x[5] <- x.cept.12
    cell.y[5] <- y.sect
    m <- (y.cept.13.2 - y.cept.13.1)/(x.cept.13.2 - x.cept.13.1)
    y0 <- y.centres[2]
    x0 <- x.centres[2]
    b <- y.cept.13.1 - m * x.cept.13.1
    x.sect <- (m * y0 + x0 - m * b)/(m^2 + 1) + sqrt(r2^2 - ((y0 -
        m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
    y.sect <- (m^2 * y0 + m * x0 + b)/(m^2 + 1) + m * sqrt(r2^2 -
        ((y0 - m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
    cell.x[3] <- (x.cept.13.1 + x.sect)/2
    cell.y[3] <- (y.cept.13.1 + y.sect)/2
    m <- (y.cept.23.2 - y.cept.23.1)/(x.cept.23.2 - x.cept.23.1)
    y0 <- y.centres[1]
    x0 <- x.centres[1]
    b <- y.cept.23.1 - m * x.cept.23.1
    x.sect <- (m * y0 + x0 - m * b)/(m^2 + 1) - sqrt(r1^2 - ((y0 -
        m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
    y.sect <- (m^2 * y0 + m * x0 + b)/(m^2 + 1) - m * sqrt(r1^2 -
        ((y0 - m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
    cell.x[1] <- (x.cept.23.1 + x.sect)/2
    cell.y[1] <- (y.cept.23.1 + y.sect)/2
    y.sect <- -sqrt(r3^2 - (x.cept.12 - x.centres[3])^2) + y.centres[3]
    cell.x[7] <- x.cept.12
    cell.y[7] <- (y.cept.12.2 + y.sect)/2
    m <- (y.cept.23.2 - y.cept.23.1)/(x.cept.23.2 - x.cept.23.1)
    y0 <- y.centres[1]
    x0 <- x.centres[1]
    b <- y.cept.23.1 - m * x.cept.23.1
    x.sect <- (m * y0 + x0 - m * b)/(m^2 + 1) + sqrt(r1^2 - ((y0 -
        m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
    y.sect <- (m^2 * y0 + m * x0 + b)/(m^2 + 1) + m * sqrt(r1^2 -
        ((y0 - m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
    cell.x[6] <- (x.cept.23.2 + x.sect)/2
    cell.y[6] <- (y.cept.23.2 + y.sect)/2
    m <- (y.cept.13.2 - y.cept.13.1)/(x.cept.13.2 - x.cept.13.1)
    y0 <- y.centres[2]
    x0 <- x.centres[2]
    b <- y.cept.13.1 - m * x.cept.13.1
    x.sect <- (m * y0 + x0 - m * b)/(m^2 + 1) - sqrt(r2^2 - ((y0 -
        m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
    y.sect <- (m^2 * y0 + m * x0 + b)/(m^2 + 1) - m * sqrt(r2^2 -
        ((y0 - m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
    cell.x[4] <- (x.cept.13.2 + x.sect)/2
    cell.y[4] <- (y.cept.13.2 + y.sect)/2
    y.sect <- sqrt(r3^2 - (x.cept.12 - x.centres[3])^2) + y.centres[3]
    cell.x[2] <- x.cept.12
    cell.y[2] <- (y.cept.12.1 + y.sect)/2
    for (i in 1:7) {
        grob.list <- gList(grob.list, textGrob(label = cell.labels[i],
            x = cell.x[i], y = cell.y[i], gp = gpar(col = label.col[i],
                cex = cex[i], fontface = fontface[i], fontfamily = fontfamily[i])))
    }
    text.location.mapping <- c(1, 3, 7)
    for (i in 1:3) {
        if ("outer" == cat.default.pos) {
            this.cat.pos <- find.cat.pos(x = x.centres[i], y = y.centres[i],
                pos = cat.pos[i], dist = cat.dist[i], r = radii[i])
        }
        else if ("text" == cat.default.pos) {
            this.cat.pos <- find.cat.pos(x = cell.x[text.location.mapping[i]],
                y = cell.y[text.location.mapping[i]], pos = cat.pos[i],
                dist = cat.dist[i])
        }
        else {
            stop("Invalid setting of cat.default.pos")
        }
        grob.list <- gList(grob.list, textGrob(label = category[i],
            x = this.cat.pos$x, y = this.cat.pos$y, just = cat.just[[i]],
            gp = gpar(col = cat.col[i], cex = cat.cex[i], fontface = cat.fontface[i],
                fontfamily = cat.fontfamily[i])))
    }
    grob.list <- VennDiagram::adjust.venn(VennDiagram::rotate.venn.degrees(gList1 = grob.list,
        angle = rotation.degree, x.centre = rotation.centre[1],
        y.centre = rotation.centre[2]), ...)
    if (ind) {
        grid.draw(grob.list)
    }
    return(grob.list)
}

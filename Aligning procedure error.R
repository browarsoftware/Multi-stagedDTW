#plot
path <- "e:\\Publikacje\\eecs 2018\\extended2_02_11\\" 
bearhug <- read.csv(paste(path, "bearhuh.txt", sep = ""), header = FALSE)
bearhug <- bearhug[1:25,]

righthandbetweenshoulders <- read.csv(paste(path, "righthandbetweenshoulders.txt", sep  = ""), header = FALSE)
righthandbetweenshoulders <- righthandbetweenshoulders[1:25,]

rightshouldertouch <- read.csv(paste(path, "rightshouldertouch.txt", sep = ""), header = FALSE)
rightshouldertouch <- rightshouldertouch[1:25,]

risingrightarmfrontally <- read.csv(paste(path, "risingrightarmfrontally.txt", sep = ""), header = FALSE)
risingrightarmfrontally <- risingrightarmfrontally[1:25,]

shoulderelevation <- read.csv(paste(path, "shoulderelevation.txt", sep = ""), header = FALSE)
shoulderelevation <- shoulderelevation[1:25,]

weightslifting <- read.csv(paste(path, "weightslifting.txt", sep = ""), header = FALSE)
weightslifting <- weightslifting[1:25,]

max_value <- max(c(bearhug, righthandbetweenshoulders, rightshouldertouch, risingrightarmfrontally, shoulderelevation,
      weightslifting))

plotlineanddotes <- function(vec, col, pch)
{
  lines(x = 1:length(vec), y = vec, 
        type = 'l', xlab = 'Iterations', ylab = 'Aligning error [cm]', col = col)
  points(x = 1:length(vec), y = vec, col = col, pch = pch)  
}

plot(x = 1:length(bearhug), y = bearhug, 
     type = 'l', xlab = 'Iterations', ylab = 'Aligning error [cm]', col = 'green', ylim = c(0, max_value * 1.1))
points(x = 1:length(bearhug), y = bearhug, col = 'green', pch=1)

plotlineanddotes(righthandbetweenshoulders, 'cyan', 2)
plotlineanddotes(rightshouldertouch, 'magenta', 3)
plotlineanddotes(risingrightarmfrontally, 'blue', 18)
plotlineanddotes(shoulderelevation, 'red', 4)
plotlineanddotes(weightslifting, 'yellow', 5)



# Add a legend
legend(15, max_value * 1.1, legend=c("Bear hug", "Hand between shoulders", "Shoulder touch", "Arm frontally",
                        "Shoulder elevation", "Weights lifting"),
       col=c("green", "cyan", "magenta", "blue", "red", "yellow"), lty=1, cex=0.8)


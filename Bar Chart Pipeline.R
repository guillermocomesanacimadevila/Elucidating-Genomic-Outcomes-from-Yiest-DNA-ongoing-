# Load ggplot2 library
library(ggplot2)

# Sample data frame
data <- data.frame(
  category = c("CUG", "UUA, UUG, CUU, CUC, CUA"),
  value = c(8244, 48409)
)

# Create a bar chart
ggplot(data, aes(x = category, y = value)) +
  geom_bar(stat = "identity") +
  labs(title = "Bar Chart of Category vs Value",
       x = "Category",
       y = "Value") +
  theme_minimal()


# Bar chart 2
data <- data.frame(
  category = c("CUG", "UUA", "UUG", "CUU", "CUC", "CUA"),
  value = c(8244, 5359, 15743, 11630, 9958, 5719)
)

# Create a bar chart
ggplot(data, aes(x = category, y = value)) +
  geom_bar(stat = "identity") +
  labs(title = "Bar Chart of Category vs Value",
       x = "Category",
       y = "Value") +
  theme_minimal()

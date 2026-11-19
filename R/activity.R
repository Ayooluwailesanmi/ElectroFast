# Single-ion activity coefficients from the Davies equation,
#   log10(gamma) = -A z^2 ( sqrt(I) / (1 + sqrt(I)) - 0.3 I ),
# with A = 0.508 at 25 C (aqion.de, "Activity Coefficients (Activity Models)",
# which gives the equation's validity as I <= 0.5 M). For standards of a 1:1 salt
# such as LiCl or LiClO4 in water with no background electrolyte (the Li+
# standards are LiClO4 in water, nothing added), the ionic strength equals
# the salt concentration.

davies_gamma = function(I, z = 1, A = 0.508) {
  10^(-A * z^2 * (sqrt(I) / (1 + sqrt(I)) - 0.3 * I))
}

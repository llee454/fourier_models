open! Core
open! Gsl

let prod2 x_range y_range f g =
  let open Range in
  (Gsl.Integrate.integration_qagp ~lower:x_range.lower ~upper:x_range.upper ~singularities:x_range.singularities ~f:(fun x ->
    (Gsl.Integrate.integration_qagp ~lower:y_range.lower ~upper:y_range.upper ~singularities:y_range.singularities ~f:(fun y ->
      Float.((f x y) * (g x y))
    )).out
  )).out

let proj2 x_range y_range f g =
  let open Float in
  (prod2 x_range y_range f g)/
  (prod2 x_range y_range g g)

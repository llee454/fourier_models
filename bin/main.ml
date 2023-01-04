open! Core
open! Aux
open! Lwt.Syntax
open! Lwt.Infix
open! Gsl

module FS_maxima = Fourier_series.Make (Fourier_series.To_maxima)

module FS_float = Fourier_series.Make (Fourier_series.To_float)

let get_trig_fn_coeffs2 n x_range y_range f =
  FS_float.get_fns2 n x_range y_range
  |> Array.map ~f:(fun g -> Orthog.proj2 x_range y_range f g)

let pdf_inv x x1 = Gsl.pdf_normal ~mean:0.0 ~std:1.0 @@ x /. x1

let main =
  Lwt_main.run @@
    let x_range  = Range.{ lower = (-1.0); upper = 1.0; singularities = [| 0.0 |] } in
    let x0_range = Range.{ lower = (-1.0); upper = 1.0; singularities = [| 0.0 |] } in
    let ks = pdf_inv |> get_trig_fn_coeffs2 3 x_range x0_range in 
    let () =
      let open FS_maxima in
      let f =
        get_fns 3 x0_range
        |> Array.mapi ~f:(fun i term -> (sprintf "as [%d]" Int.(i + 1)) * term "x0")
        |> sum
      in
      let g =
        get_fns2 3 x_range x0_range
        |> Array.mapi ~f:(fun i term -> (const ks.(i)) * term "x" "x0")
        |> sum
      in
      f * g |> printf {maxima|
        trigexpand: true$
        ratprint: false$
        f (x) := trigreduce (expand (integrate (trigexpand (%s), x0, -1, 1)))$
      |maxima}
    in
    Lwt.return_unit

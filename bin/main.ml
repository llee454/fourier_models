open! Core
open! Lwt.Syntax
open! Lwt.Infix
open! Gsl

module Range = struct
  type t = {
    lower: float;
    upper: float;
    singularities: float array
  } [@@deriving fields]
end

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

(**
  Returns an array of trigonemtric functions for use in a truncated
  fourier series. Accepts an argument, n, which determines the
  maximum degree for the series. The functions are returned as 1
  followed by a sequence of "blocks." Each block has a prefix of
  the form: [cos y; sin y; cos x; cos y] followed by a "body." Each
  body consists of a sequence of "squares" of the form: [
  cos (a x) cos (b y); cos (a x) sin (b y); sin (a x) cos (b y);
  sin (a x) sin (b y)] where a and b are determined by [(a, b)] =
  [(1, n); (n, 1); (2, n); (n, 2); ...] ending in (n, n).

  So for example, when n = 2 we get:
  [1;
   cos x; sin x; cos y; sin y; # block 1 prefix
   <nothing> # block 1 squares 
   cos x cos y; cos x sin y; sin x cos y; sin x sin y; # block 1 suffix
   cos 2x; sin 2x; cos 2y; sin 2y; # block 2 prefix
   cos x cos 2y; cos x sin 2y; sin x cos 2y; sin x sin 2y; # block 2 square 1
   cos 2x cos y; cos 2x sin y; sin 2x cos y; sin 2x sin y;
   cos 2x cos 2y; cos 2x sin 2y; sin 2x cos 2y; sin 2x sin 2y # block 2 suffix
  ] 
*)
let get_trig_fns n x_range y_range =
  let open Range in
  let open Float in
  let x_delta, y_delta = x_range.upper - x_range.lower, y_range.upper - y_range.lower in
  let fcos delta k x = cos (float k * 2.0 * pi * x / delta) in
  let fsin delta k x = sin (float k * 2.0 * pi * x / delta) in
  let xs = Queue.create ~capacity:Int.(1 + 4*n + 4*(n*n)) () in
  Queue.enqueue xs (fun _ _ -> 1.0);
  for block_idx = 1 to n do
    Queue.enqueue_all xs [
      (fun x _y -> fcos x_delta block_idx x);
      (fun x _y -> fsin x_delta block_idx x);
      (fun _x y -> fcos y_delta block_idx y);
      (fun _x y -> fsin y_delta block_idx y)
    ];
    for i = 1 to Int.(block_idx - 1) do
      Queue.enqueue_all xs [
        (fun x y -> fcos x_delta i x * fcos y_delta block_idx y);
        (fun x y -> fcos x_delta i x * fsin y_delta block_idx y);
        (fun x y -> fsin x_delta i x * fcos y_delta block_idx y);
        (fun x y -> fsin x_delta i x * fsin y_delta block_idx y);
        (fun x y -> fcos x_delta block_idx x * fcos y_delta i y);
        (fun x y -> fcos x_delta block_idx x * fsin y_delta i y);
        (fun x y -> fsin x_delta block_idx x * fcos y_delta i y);
        (fun x y -> fsin x_delta block_idx x * fsin y_delta i y);
      ]
    done;
    Queue.enqueue_all xs [
      (fun x y -> fcos x_delta block_idx x * fcos y_delta block_idx y);
      (fun x y -> fcos x_delta block_idx x * fsin y_delta block_idx y);
      (fun x y -> fsin x_delta block_idx x * fcos y_delta block_idx y);
      (fun x y -> fsin x_delta block_idx x * fsin y_delta block_idx y)
    ];
  done;
  Queue.to_array xs

let get_trig_fn_coeffs n x_range y_range f =
  get_trig_fns n x_range y_range
  |> Array.map ~f:(fun g -> proj2 x_range y_range f g)

let pdf_inv x x1 =
  let open Float in
  Gsl.pdf_normal ~mean:0.0 ~std:1.0 @@ x / x1 

let main =
  Lwt_main.run @@
    let* () =
      let ks = pdf_inv
      |> get_trig_fn_coeffs 3
         { lower = (-1.0); upper = 1.0; singularities = [| 0.0 |] }
         { lower = (-1.0); upper = 1.0; singularities = [| 0.0 |] }
      in 
      printf "ks: [\n";
      for i = 0 to Array.length ks - 1 do
        printf !"%0.08f, " ks.(i)
      done;
      printf "]$\n";
      Lwt.return_unit
    in
(*
    let Linear_fit.{c0; c1} = Linear_fit.f [|0.0; 1.0; 2.0|] [|5.0; 10.0; 20.0|] in 
    let* () = Lwt_io.printlf "linear regression: %f + %f x" c0 c1 in
    let Integrate.{ out; err; _ } = Integrate.f ~f:(fun x -> Float.log1p x) ~lower:0.1 ~upper:0.5 in
    let* () = Lwt_io.printlf "integration: result %f err %f" out err in
    let result =
      Simulated_annealing.f {
        copy = (fun x -> ref !x);
        energy = (fun x -> (!x +. 7.0) *. (!x +. 7.0));
        step = (fun x dist -> x := !x +. dist);
        dist = (fun x y -> Float.abs (!x -. !y));
        init = (ref 5.0);
        print = None
      }
    in
    let* () = Lwt_io.printlf "simulated annealing result: %f" !result in
*)
    Lwt.return_unit

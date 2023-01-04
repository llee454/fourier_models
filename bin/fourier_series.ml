open! Core
open! Gsl

module type Intf = sig
  type t

  val ( + ) : t -> t -> t

  val ( * ) : t -> t -> t

  val ( / ) : t -> t -> t

  val pi : t

  val sin : t -> t

  val cos : t -> t

  val sum : t array -> t

  val const : float -> t
end

module To_maxima : Intf with type t = string = struct
  type t = string

  let ( + ) = sprintf "(%s + %s)"

  let ( * ) = sprintf "(%s * %s)"

  let ( / ) = sprintf "(%s / %s)" 

  let pi = "%pi"

  let sin = sprintf "sin (%s)"

  let cos = sprintf "cos (%s)"

  let sum = Array.fold ~init:"" ~f:(fun expr term -> expr + term)

  let const = sprintf "%f"
end

module To_float : Intf with type t = float = struct
  include Float

  let const = Fn.id

  let sum = Array.fold ~init:0.0 ~f:(fun expr term -> expr + term)
end

module Make (M : Intf) = struct
  include M

  (**
    Returns the first n-th order terms terms in a truncated fourier
    series scaled to approximate a function over the given range.
  *)
  let get_fns n (range : Range.t) =
    let open M in
    let delta = Float.(range.upper - range.lower) in
    let fcos k x = cos (const (float k) * const 2.0 * pi * x / const delta) in
    let fsin k x = sin (const (float k) * const 2.0 * pi * x / const delta) in
    let xs = Queue.create ~capacity:Int.(2 * n) () in
    Queue.enqueue xs (Fn.const (const 1.0));
    for i = 1 to n do
      Queue.enqueue_all xs [fcos i; fsin i]
    done;
    Queue.to_array xs

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
  let get_fns2 n x_range y_range =
    let open Range in
    let open M in
    let x_delta, y_delta = Float.(x_range.upper - x_range.lower, y_range.upper - y_range.lower) in
    let fcos delta k x = cos (const (float k) * const 2.0 * pi * x / const delta) in
    let fsin delta k x = sin (const (float k) * const 2.0 * pi * x / const delta) in
    let xs = Queue.create ~capacity:Int.(1 + 4*n + 4*(n*n)) () in
    Queue.enqueue xs (fun _ _ -> const 1.0);
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
end

open! Core

type t = {
  lower: float;
  upper: float;
  singularities: float array
} [@@deriving fields]

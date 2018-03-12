(********** POLYNÔMIAL ***********)

module type IPoly =
	sig
		
		type poly
		
  	val poly_empty : poly
  	val poly_unit : poly

    val poly_n : int -> poly
		val poly_max_degree : poly -> int
		val poly_first_degree : poly -> int
		val poly_tl : poly -> poly
		val poly_contains : poly -> int -> bool

    val add_poly : poly -> poly -> poly
		val mult_poly : poly -> poly -> poly
    val mult_karatsuba_poly : poly -> poly -> poly
    val newton_poly : poly -> int -> poly
    val quomod_poly : poly -> poly -> poly * poly
    val quomod_newton_poly : poly -> poly -> poly * poly
    val quomod_asc_poly : poly -> poly -> int -> poly * poly
    val euclide_poly : poly -> poly -> poly
    val order_poly : poly -> int

    val irreducible_poly : int -> poly
    val primitive_poly : int -> poly
		val poly_random_not_empty : int -> poly

    val poly_of_string : string -> poly
    val string_of_poly : poly -> string
		
	end

module Poly : IPoly =
	struct

		(* ---------- Mathematical --------- *)

    (* Renvoie le logarithme à base 2 (arrondie au supérieur) *)
    let log2 = function n ->
			int_of_float (ceil (log (float_of_int n) /. log 2.))
		;;

		(* Renvoie a à la puissance b *)
		let pow = fun a b ->
			int_of_float ((float_of_int b) ** (float_of_int (a)))
		;;

    (* ----------            ---------- *)

		(* ---------- Polynomial ---------- *)

		type poly = int list;;

		let poly_empty = [];;

		let poly_unit  = [0];;

		(* Renvoie le monôme de degré n *)
		let poly_n = function n ->
			[n]
		;;

		(* Renvoie le succession immédiat du polynôme p *)
		let poly_succ = function p ->
			let rec aux = fun p n ->
				match p with
				| [] -> [n]
				| d::pp when d = n -> (aux pp (n + 1))
				| d::pp -> (n::p)
			in (aux p 0)
		;;

		(* Renvoie le degré maximal du polynôme
   		 	si le polynôme est dans le bon ordre
    		sinon renvoie -1 *)
    let poly_max_degree = function a ->
			if a = poly_empty then -1 else List.hd (List.rev a)
		;;

    (* Renvoie le premier degrée de la liste
    		si la liste est vide renvoie -1 *)
    let poly_first_degree = function a ->
    	match a with
    	| d::_ -> d
    	| _ -> -1
    ;;

		(* Renvoie tout le polynôme sauf le degré minimal non nul du polynôme p *)
		let poly_tl = function p ->
			match p with
			| _::ll -> ll
			| _ -> failwith "poly_tl: poly_empty"
		;;

    (* Coupe le polynôme a à partir du degrée n *)
    let split_at_poly = fun a n ->
    	let rec aux = fun a l h ->
    		match a with
    		| [] -> (List.rev l, List.rev h)
    		| d::ll when d <= n -> (aux ll (d::l) h)
    		| d::ll -> (aux ll l ((d - n)::h))
    	in (aux a poly_empty poly_empty)
    ;;

    (* Renvoie le renversé d'ordre n du polynôme p *)
    let renv_poly = fun p n ->
    	if n < (poly_max_degree p) then	failwith "renv_poly: wrong degree"
    	else let rec aux = fun p acc ->
    		match p with
  			| [] -> acc
  			| d::pp -> (aux pp (n - d::acc))
  		in (aux p poly_empty)
    ;;

		(* Effectue le module n sur le polynôme p *)
    let mod_poly = fun p n ->
    	let rec aux = fun p acc ->
    		match p with
    		| d::pp when d < n -> (aux pp (d::acc))
    		| _ -> acc
    	in (aux p poly_empty)
    ;;

    (* Renvoie la moitié du nombre et rajoute 1 si celui-ci est impair *)
    let half_degree = function n ->
    	let n2 = (n / 2) in
    		if (n2 mod 2) = 0 then n2
    		else (n2 + 1)
    ;;

		(* Teste si le polynôme contient un terme de degré n *)
		let poly_contains = fun p n ->
			let rec aux = function p ->
				match p with
				| d::_ when (d = n) -> true
				| d::pp when (d < n) -> (aux pp)
				| _ -> false
			in (aux p)
		;;

		(* ----------            ---------- *)

		(* ---------- Operations ---------- *)

    (* Effectue l'addition entre deux polynômes p1 et p2 *)
    let add_poly = fun p1 p2 ->
    	let rec aux = fun p1 p2 acc ->
      	match (p1, p2) with
      	| [], [] -> List.rev acc
    		| d::pp, [] -> (aux pp [] (d::acc))
    		| [], d::pp -> (aux [] pp (d::acc))
    		| d1::pp1, d2::pp2 when d1 > d2 -> (aux p1 pp2 (d2::acc))
    		| d1::pp1, d2::pp2 when d1 < d2 -> (aux pp1 p2 (d1::acc))
    		| d1::pp1, d2::pp2 -> (aux pp1 pp2 acc)
    	in (aux p1 p2 [])
    ;;

		(* Renvoie un polynome au hasard de degré plus petit où égal à n et non nul *)
		let poly_random_not_empty = function n ->
			let rec aux = fun k acc ->
				match k with
				| k when (k >= (n + 1)) -> acc
				| k -> let r = (Random.int 2)
					in if (r = 0) then (aux (k + 1) acc)
					else (aux (k + 1) (add_poly (poly_n k) acc))
			in (aux 1 poly_unit)
		;;

    (* Effectue la multiplication naive entre deux polynômes p1 et p2 *)
    let mult_poly = fun p1 p2 ->
    	let rec aux = fun a b acc ->
    		match (a, b) with
    		| ([], []) | ([], _) -> acc
    		| _::pp, [] -> (aux pp p2 acc)
    		| d1::_, d2::pp -> (aux a pp (add_poly acc [d1 + d2]))
    	in (aux p1 p2 [])
    ;;

		(* Ajoute n à tous les degrés du polynôme p *)
    let mult_xn_poly = fun n p ->
    	let rec aux = fun p acc ->
      	match p with
      	| [] -> acc
      	| d::pp -> aux pp ((d + n)::acc)
    	in List.rev (aux p [])
    ;;

		(* Effectue la multiplication de p1 et p2 par la méthode de Karatsuba *)
    let rec mult_karatsuba_poly = fun p1 p2 ->
    	if (max (poly_max_degree p1) (poly_max_degree p2)) < 3 then
    		(mult_poly p1 p2)
    	else let m2 = (half_degree (max (poly_max_degree p1) (poly_max_degree p2)))
    		in let (l1, h1) = (split_at_poly p1 m2)	and (l2, h2) = (split_at_poly p2 m2)
    		in let	z0 = (mult_karatsuba_poly l1 l2)
    		and			z1 = (mult_karatsuba_poly (add_poly l1 h1) (add_poly l2 h2))
    		and			z2 = (mult_karatsuba_poly h1 h2)
    		in (add_poly
    				(add_poly (mult_xn_poly (2 * m2) z2)
    									(mult_xn_poly m2 (add_poly (add_poly z1 z2) z0)))
    				z0)
    ;;

		(* Calcule le newton n du polynôme p *)
    let newton_poly = fun p n ->
    	if (poly_first_degree p != 0) then failwith "newton_poly: the constant term is not 1"
    	else let rec aux = fun g r i ->
    		if (i > r) then g
				else (aux (mult_karatsuba_poly p (mult_karatsuba_poly g g)) r (i + 1))
  		in (aux poly_unit (log2 n) 1)
    ;;

		let quomod_poly = fun p1 p2 ->
			let rec aux = fun q r ->
				if (r = poly_empty || (poly_max_degree r) < (poly_max_degree p2)) then	(q, r)
				else let t = (poly_n ((poly_max_degree r) - (poly_max_degree p2)))
					in (aux (add_poly q t) (add_poly r (mult_karatsuba_poly t p2)))
			in (aux poly_empty p1)
		;;

		(* Effectue la division de deux polynômes p1 et p2 *)
		(* 	renvoie le quotient et le reste *)
    let quomod_newton_poly = fun p1 p2 ->
    	let n = (poly_max_degree p1) and m = (poly_max_degree p2)
    	in if (n < m) then (poly_empty, p1)
    	else let p = (newton_poly (renv_poly p2 m) (n - m + 1))
    		in let q = (renv_poly (mod_poly (mult_karatsuba_poly (renv_poly p1 n) p) (n - m + 1)) (n - m))
    		in (List.rev q, (add_poly p1 (mult_karatsuba_poly p2 q)))
    ;;

		(* Effectue la division des puissances croissantes de deux polynômes p1 et p2 *)
		let quomod_asc_poly = fun p1 p2 n ->
			let rec aux = fun p1 acc ->
  			if ((poly_first_degree p1) > n) then (acc, p1)
				else let q = (poly_first_degree p1) - (poly_first_degree p2)
					in (aux (add_poly (mult_xn_poly q p2) p1) (add_poly (poly_n q) acc))
			in (aux p1 poly_empty)
		;;

		(* Donne le PGCD de deux polynômes p1 et p2 *)
    let rec euclide_poly = fun p1 p2 ->
    	let (_, r) = (quomod_poly p1 p2)
			in if (r = poly_empty) then p2
				else (euclide_poly p2 r)
    ;;

		(* Renvoie l'ordre du polynôme p *)
		let order_poly = function p ->
    	if poly_max_degree p <= 1 then failwith "order_poly: degree <= 1"
    	else if (poly_first_degree p != 0) then failwith "order_poly: the constant term is not 1"
    		else let rec aux = function n ->
      		let (_, r) = (quomod_poly (poly_n n) p)
      		in if (r = poly_unit) then n
      			else (aux (n + 1))
      	in (aux 1)
  	;;

		(* Test pour savoir si le polynôme p est irréductible *)
		let is_irreducible_poly = function p ->
			let rec aux = function acc ->
				let (_, r) = (quomod_poly p acc) and succ = (poly_succ acc)
				in if (r != poly_empty) then
					if (poly_max_degree succ = (poly_max_degree p)) then true
					else (aux succ)
				else false
			in (aux (poly_n 1))
		;;

		(* Renvoie le premier polynôme irréductible de degré n *)
		let irreducible_poly = function n ->
			if (n <= 0) then failwith "irreducible_poly: n <= 0"
			else let rec aux = function acc ->
				if (is_irreducible_poly acc) then acc
				else (aux (poly_succ acc))
			in (aux (poly_n n))
		;;

		(* Test pour savoir si le polynôme p est primitif *)
		let is_primitive_poly = function p ->
			let n = ((pow (poly_max_degree p) 2) - 1)
			in let (_, r) = (quomod_poly (poly_n n) p)
			in ((r = poly_unit) && (is_irreducible_poly p) && (order_poly p = n))
		;;

		(* Renvoie le premier polynôme primitif de degré n *)
		let primitive_poly = function n ->
			if n <= 0 then failwith "primitive_poly: n <= 0"
			else let rec aux = function acc ->
				if (is_primitive_poly acc) then acc
				else (aux (poly_succ (poly_succ acc)))
			in (aux (poly_succ (poly_n n)))
		;;

		(* ----------           ---------- *)

		(* ---------- Tools ---------- *)

		let char_unknown = 'X';;

		let char_numbers = ['0'; '1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'];;

		let good_char = char_unknown::char_numbers;;

		(* Obtient le reste de la chaine s en enlevant le premier caractère *)
		let tl_string = function s ->
			(String.sub s 1 (String.length s - 1))
		;;

		(* Renvoie la chaîne de caractère constitué
				de tout les chiffres successifs dans s
				renvoie aussi le reste de cette chaine *)
		let number_string = function s ->
			let rec aux = fun s snum ->
				if (String.length s = 0 || (not (List.mem s.[0] char_numbers))) then
					if (snum = "") then failwith "number_string: bad format"
					else (s, int_of_string snum)
				else (aux (tl_string s) (snum^(Char.escaped s.[0])))
			in (aux s "")
		;;

		(* Nettoie la chaîne de caractère :
		 		supprime les caractères non correcte
				seul le caractère X est accepté mais aussi les chiffres *)
		let clean_string = function s ->
			let rec aux = fun s acc ->
				if String.length s = 0 then	acc
				else let c  = s.[0] and ss = (String.sub s 1 (String.length s - 1))
					in if (List.mem c good_char) then	(aux ss (acc^(Char.escaped c)))
					else (aux ss acc)
			in (aux (String.uppercase s) "")
		;;

		(* ----------       ---------- *)

		(* ---------- Conversion ---------- *)

		let poly_of_string = function s ->
			let rec aux = fun s acc ->
				if String.length s = 0 then	(List.sort_uniq (fun a b -> a - b) acc)
				else let c = s.[0]
					in if (c = char_unknown) then
						let (ss, d) = number_string (tl_string s)
						in aux ss (d::acc)
					else if (c = '1') then (aux (tl_string s) (0::acc))
					else failwith "poly_of_string: bad format"
			in (aux (clean_string s) [])
		;;

		let string_of_poly = function p ->
			let rec aux = fun p acc ->
				match p with
				| [] -> acc
				| d::[] -> "X^"^(string_of_int d)^acc
				| d::ll -> (aux ll (" + X^"^(string_of_int d)^acc))
			in (aux p "")
		;;

		(* ----------            ---------- *)
	end


(********** LFRS ***********)

open Poly

module type ILFSR =
	sig
		
    type lfsr

    val lf1 : lfsr
    val lf2 : lfsr

    val rn_list_lfsr : lfsr -> int -> int list
    val rn_lfsr : lfsr -> int -> int

		val poly_of_lfsr : lfsr -> int * Poly.poly * Poly.poly
		val lfsr_of_poly : int * Poly.poly * Poly.poly -> lfsr
    val min_lfsr : lfsr -> lfsr
		val good_lfsr : int -> lfsr

    val string_of_lfsr : lfsr -> string
	
	end

module LFSR : ILFSR =
	struct

		(* ---------- Mathematical --------- *)

		let pow = fun a b ->
			int_of_float ((float_of_int b) ** (float_of_int (a)))
		;;

    (* ----------            ---------- *)

		type lfsr = {registre : poly; l : int; branchement : poly};;

		(* ----------    Exemples  ---------- *)

		let lf1 = {registre = (poly_of_string "1 + x^3 + x^6 + x^9");
											l = 10 ;
						branchement = (poly_of_string "x^1 + x^3 + x^4 + x^7 + x^10")};;

    let lf2 = {registre = (poly_of_string "1");
											l = 3 ;
						branchement = (poly_of_string "x^3")};;

    (* ----------            ---------- *)

    (* ---------- Opérations ---------- *)

		(* Renvoie l'élément (n - 1) de la liste l *)
		let elem_n_in_list = fun n l ->
			let rec aux = fun l acc ->
				match l with
				| e::ll when (acc = n) -> e
				| e::ll -> (aux ll (acc + 1))
				| _ -> failwith "elem_n_in_list : list is too tiny."
			in (aux l 0)
		;;

		(* Renvoie le R(n + 1) selon le lfsr l *)
		let lfsr_add_next = fun l l_pred ->
			let rec aux = fun p acc ->
				match p with
				| p when (p = poly_empty) -> (List.rev ((acc mod 2)::(List.rev l_pred)))
				| p -> (aux (poly_tl p)
										((elem_n_in_list ((List.length l_pred) - (poly_first_degree p)) l_pred) + acc))
			in (aux l.branchement 0)
		;;

		(* Renvoie le coeff du degré n du registre du lfsr l *)
		let registre_lfsr = fun l n ->
			if (poly_contains l.registre n) then 1 else 0
		;;

		(* Renvoie une liste des premières valeurs selon le registre du lfsr l *)
		let ini_lfsr = fun l n ->
			let rec aux = fun k acc ->
				match k with
				| k when (k = (n + 1)) -> (List.rev acc)
				| k -> (aux (k + 1) ((registre_lfsr l k)::acc))
			in (aux 0 [])
		;;

		(* Renvoie une liste des n des premières valeurs du lfsr l *)
		let rn_list_lfsr = fun l n ->
			if (n < l.l) then (ini_lfsr l n)
			else let ini = (ini_lfsr l (l.l - 1))
			in let rec aux = fun k acc ->
				if (k = (n + 1)) then acc
				else (aux (k + 1) (lfsr_add_next l acc))
			in (aux l.l ini)
		;;

		(* Renvoie la valeur n du lfsr l *)
		let rn_lfsr = fun l n ->
			(List.hd (List.rev (rn_list_lfsr l n)))
		;;

		(* Renvoie le coeff du degré n du polynôme de rétroaction du lfsr l *)
		let rx_lfsr = fun l n ->
			if (poly_contains (add_poly l.branchement poly_unit) n) then 1 else 0
		;;

		(* Renvoie la somme des α(i−j) * r(j) *)
		let sum_for_g = fun l i ->
			let rec aux = fun j acc ->
				match j with
				| j when j = (i + 1) -> (acc mod 2)
				| j -> (aux (j + 1) (((rx_lfsr l (i - j)) * (rn_lfsr l j)) + acc))
			in (aux 0 0)
		;;

		(* Renvoie g(x) du lfsr l *)
		let calcul_g = function l ->
			let rec aux = fun i acc ->
				match i with
				| i when i = (l.l) -> acc
				| i -> let new_acc = if ((sum_for_g l i) = 0) then acc
						else (add_poly (poly_n i) acc)
					in (aux (i + 1) new_acc)
			in (aux 0 poly_empty)
		;;

		(* Renvoie un triplet (longueur, g(x), rétroaction) du lfsr l *)
		let poly_of_lfsr = function l ->
			let r = (add_poly l.branchement poly_unit) in
			((poly_max_degree r), (calcul_g l), r)


		(* Renvoie un lfsr à partir du triplet (l, g, r) *)
		let lfsr_of_poly = function (l, g, r) ->
			let registre = (fst (quomod_asc_poly g r (l - 1))) in
			{registre = registre; l = l; branchement = (add_poly r poly_unit)}
		;;

		(* Renvoie un triplet réduit à partir du triplet (l, g, r) *)
		let less_poly_idem_value = function (l, g, r) ->
			let t = (euclide_poly g r)
			in let new_g = (fst (quomod_newton_poly g t))
			and new_r = (fst (quomod_newton_poly r t))
			in ((poly_max_degree new_r), new_g, new_r)
		;;

    (* Renvoie un lfsr réduit à partir du lfsr l de même flux de valeurs *)
		let min_lfsr = function l ->
			(lfsr_of_poly (less_poly_idem_value (poly_of_lfsr l)))
		;;

		(* Renvoie un 'bon' lfsr *)
		let good_lfsr = function n ->
			if (n < 1) then failwith "good_lfsr : n < 1"
			else let r = (primitive_poly n)
			and g = (poly_random_not_empty (n - 1))
			in (lfsr_of_poly (n, g, r))
		;;

		(* ----------            ---------- *)

		(* ---------- Conversion ---------- *)

		(* Retourne le lfsr sous forme de chaîne de caractère *)
		let string_of_lfsr = function l ->
			"LFSR : l = "^ (string_of_int l.l)
			^ "; *" ^ (string_of_poly l.registre;)
			^ "*; *" ^ (string_of_poly l.branchement)
			^ "*"
		;;

		(* ----------            ---------- *)
	end

open LFSR;;


module type ICrypto =
	sig
		
		type encrypted
		
		val encrypt : LFSR.lfsr -> string -> encrypted
    val decrypt : LFSR.lfsr -> encrypted -> string
		val string_of_encrypted : encrypted -> string
		
	end

module Crypto : ICrypto =
	struct
		
		type encrypted = int list;;

		(* Renvoie le xor entre r(k) du lfsr et b *)
		let xor_lfsr = fun lfsr k b ->
			(((rn_lfsr lfsr k) + b) mod 2)
		;;

		(* Renvoie les xor entre les r(k) du lfsr et les élèments de l *)
		let xor_bits_lfsr = fun lfsr l ->
			let rec aux = fun l k acc ->
				match l with
				| [] -> (List.rev acc)
				| b::ll -> (aux ll (k + 1) ((xor_lfsr lfsr k b)::acc))
			in (aux l 0 [])
		;;


		(* Convertie un entier en binaire sur b bits *)
		let bits_of_int = fun n base ->
			let rec aux = fun k acc ->
				match k with
				| 0 when ((List.length acc) >= base) -> acc
				| 0 -> (aux 0 (0::acc))
				| k -> (aux (k/2) ((k mod 2)::acc))
			in (aux n [])
		;;

		(* Convertie un entier en binaire *)
		let int_of_bits = function l ->
			let rec aux = fun l acc ->
				match l with
				| [] -> acc
				| k::ll -> (aux ll (acc * 2 + k))
			in (aux l 0)
		;;

		(* Convertie un caractère en binaire de base bits *)
		let bits_of_char c base = (bits_of_int (int_of_char c) base);;

		(* Convertie un nombre binaire en caractère *)
		let char_of_bits b = (char_of_int (int_of_bits b));;

		(* Encode un caractère en un nombre binaire *)
		let encrypt_char = fun lfsr c ->
			(xor_bits_lfsr lfsr (bits_of_char c 7))
		;;

		(* Décode un nombre binaire en un caractère *)
		let decrypt_bits = fun lfsr b ->
			(char_of_bits (xor_bits_lfsr lfsr b))
		;;

		(* Donne la liste du caractère de tête à n et la liste restante *)
		let list_cut_head = fun l n ->
			if ((List.length l) < n) then failwith "list_cut_head : n > length of l"
			else let rec aux = fun l n acc ->
				match n with
				| 0 -> ((List.rev acc), l)
				| _ -> (aux (List.tl l) (n - 1) ((List.hd l)::acc))
			in (aux l n [])
		;;

		(* Encode une chaîne de caractère en un message encodé *)
		let encrypt = fun lfsr s ->
			let rec aux = fun s acc ->
				match s with
				| "" -> acc
				| s -> (aux (String.sub s 1 ((String.length s) - 1))
										(acc@(encrypt_char lfsr (s.[0]))))
			in (aux s [])
		;;

		(* Décode un message encodé en une chaîne de caractère *)
		let decrypt = fun lfsr b ->
			let rec aux = fun b acc ->
				match b with
				| [] -> acc
				| _ -> let (c, ll) = (list_cut_head b 7)
			 			in (aux ll (acc ^ (Char.escaped (decrypt_bits lfsr c))))
			in (aux b "")
		;;

		(* Renvoie le nombre binaire encodée sous forme d'une chaîne de caractère *)
		let string_of_encrypted = function b ->
				let rec aux = fun b acc ->
					match b with
					| [] -> acc
					| _ -> let (c, ll) = (list_cut_head b 7)
				 			in (aux ll (acc ^ (Char.escaped (char_of_bits c))))
				in (aux b "")
			;;

	end

open Crypto

(* Transforme un liste d'entier en une représentation
de liste en chaîne de caractère *)
let string_of_int_list = function l ->
	let rec aux = fun l acc ->
  	match l with
		| [] -> "[]"
		| e::[] -> (acc^(string_of_int e)^"]")
		| e::ll -> (aux ll (acc^(string_of_int e)^"; "))
	in (aux l "[")

(********** Test pour Poly **********)

let p1 = poly_of_string "1 + x^2 + x^5";;
let p2 = poly_of_string "X^1 + x^2 + x^4";;

Printf.printf "\n";
Printf.printf "********* Poly *********\n";
Printf.printf "p1 = %s\n" (string_of_poly p1);
Printf.printf "p2 = %s\n" (string_of_poly p2);
Printf.printf "p1 + p2 = %s\n" (string_of_poly (add_poly p1 p2));
Printf.printf "p1 x p2 = %s (naïve)\n" (string_of_poly (mult_poly p1 p2));
Printf.printf "p1 x p2 = %s (Karatsuba)\n" (string_of_poly (mult_karatsuba_poly p1 p2));
Printf.printf "p1 / p2 = %s (naïve)\n" (string_of_poly (fst (quomod_poly p1 p2)));
Printf.printf "p1 / p2 = %s (Newton)\n" (string_of_poly (fst (quomod_newton_poly p1 p2)));
Printf.printf "p1 %% p2 = %s (naïve)\n" (string_of_poly (snd (quomod_poly p1 p2)));
Printf.printf "p1 %% p2 = %s (Newton)\n" (string_of_poly (snd (quomod_newton_poly p1 p2)));
Printf.printf "p1.pgcd(p2) = %s\n" (string_of_poly (euclide_poly p1 p2));
Printf.printf "*********      *********\n";;

(**********                **********)

(********** Test pour LFSR **********)

let (l,g,r) = (poly_of_lfsr lf1);;
let (_,g2,r2) = (poly_of_lfsr lf2);;

Printf.printf "\n";
Printf.printf "********* LFSR *********\n";
Printf.printf "lf1 = %s\n" (string_of_lfsr lf1);
Printf.printf "lf2 = %s\n" (string_of_lfsr lf2);
Printf.printf "20 premières valeurs de lf1 = %s\n" (string_of_int_list (rn_list_lfsr lf1 20));
Printf.printf "20 premières valeurs de lf2 = %s\n" (string_of_int_list (rn_list_lfsr lf2 20));
Printf.printf "20-ième valeur de lf2 = %d\n" (rn_lfsr lf2 20);
Printf.printf "G(x) de lf1 = %s\n" (string_of_poly g);
Printf.printf "R(x) de lf1 = %s\n" (string_of_poly r);
Printf.printf "G(x) de lf1 = %s\n" (string_of_poly g2);
Printf.printf "R(x) de lf1 = %s\n" (string_of_poly r2);
Printf.printf "lf1 depuis G(x) et R(x) = %s\n" (string_of_lfsr (lfsr_of_poly (l, g, r)));
Printf.printf "lf1 réduit donnant lf2 = %s\n" (string_of_lfsr (min_lfsr lf1));
Printf.printf "un 'bon' lfsr aléatoire de degré 5 = %s\n" (string_of_lfsr (good_lfsr 5));
Printf.printf "*********      *********\n";;

(**********                **********)

(********** Test pour Crypto **********)

let glf = (good_lfsr 6);;
let message = "Si pour vous, votre vie ne vaut pas plus que celles des autres;"
	^ " signez votre carte de donneur maintenant et tuez vous ensuite.";;
let encrypted = (encrypt glf message);;
Printf.printf "\n";
Printf.printf "********* Crypto *********\n";
Printf.printf "glf est un 'bon' lfsr\n";
Printf.printf "Message clair     : %s\n" message;
Printf.printf "Message chiffré   : %s\n" (string_of_encrypted encrypted);
Printf.printf "Message déchiffré : %s\n" (decrypt glf encrypted);
Printf.printf "*********        *********\n";;

(**********                  **********)
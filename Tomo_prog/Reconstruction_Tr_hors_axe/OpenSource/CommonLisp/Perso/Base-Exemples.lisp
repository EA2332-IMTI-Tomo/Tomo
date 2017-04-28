;;GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
;;Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net)


;;;-----------------------------------------------------------------------------
;;; *Autodoc
;;;-----------------------------------------------------------------------------

;; search for function declaration from part of his name
(apropos "ensym")
(apropos "ensym" :cl)

;; display info about a standard or defined function
(describe 'mafonction)

;;only the "doc" field
(documentation 'mapcan 'function)

(inspect 'mafonction)


;;;-----------------------------------------------------------------------------
;;; *Basique
;;;-----------------------------------------------------------------------------

;;** CF section *Listes


(identity 'rty)
;;rty

(apply #'+ '(1 2 3))
;;6

;; idem mais ne requiert pas arguments dans une liste.
(funcall #'+ 1 2 3)
;;6

(function +)
;;#<Function + {1022B801}>
(lambda (x) (+ 100 x))
;;#<Interpreted Function (LAMBDA (X) (+ 100 X)) {810F9E71}>
((lambda (x) (+ 100 x)) 1)
;;101
(funcall #'(lambda (x) (+ 100 x)) 1)


(mkstr 'toto 3 " 5")
;;"TOTO3 5"


(defun use-lambda (a b myfun)
  (funcall myfun a b))
;;(use-lambda 1 2 #'+)
;;3
;;(use-lambda 1 2 (lambda (a b) (+ a b)))
;;(use-lambda 1 2 #'(lambda (a b) (+ a b)))
;;3


;;;-----------------------------------------------------------------------------
;;; * Input/Output
;;;-----------------------------------------------------------------------------

;;keywds: interaction keyboard input output display


;;display (t: stdout)
;;can be used for strings
(format t "~2,4F ~% ~D ~A ~2,'0D" (/ 2 3) 5 "toto" 5)
;;.6667 
;; 5 toto 05
;;NIL
;; The format string can be replaced by a string object

;;other examples: http://dept-info.labri.u-bordeaux.fr/~strandh/Teaching/Programmation-Symbolique/Common/HyperSpec/Body/sec_22-3-11.html


;; return a string (and not nil)
(format nil "toto is ~W" 4)
;;"toto is 4"

;; format strings evolved
(format nil "~3,'0D" 1)

;;from hyperspec
(defun sqrt-advisor ()
   (loop (format t "~&Number: ")
         (let ((n (parse-integer (read-line) :junk-allowed t)))
           (when (not n) (return))
           (format t "~&The square root of ~D is ~D.~%" n (sqrt n)))))



;;;-----------------------------------------------------------------------------
;;; *Comparaisons
;;;-----------------------------------------------------------------------------
 

;;eq
;;symboles identiques: teste l'adresse memoire (le reader affecte la meme adresse aux memes symboles) 

;;eql
;; aussi: nombres de meme type, caracteres

;;equal:
;; objets complexes decortiques. chaines: case sensitive
;; strings, listes, etc?

;; equalp
;; arrays, structures, hash-tables
;; elements compares un a un, case insensitive, type insensitive pour les nombres.

;;voir aussi: fonctions specialisees (rapides) typees
;;char= char/= char< char<= char-equal char-not-equal        
;;string...
;;tree-equal
;;nombres: = /= < > <= ...


;;;-----------------------------------------------------------------------------
;;; *Iterateurs standard
;;;-----------------------------------------------------------------------------


;;*1: DO

;;do fonctionne comme let: les variables declarées le sont dans un envt local
;;ce qui n'empeche pas d'utiliser des variables du toplevel.

;;exemple basique
(do ((i 1 (1+ i))
     (j 0 (1- j))
     (k 0))
    ((> (- i j) 5) (list i j k))
  (progn (setf k (+ 1 k))
	 (format t "~% ~A ~A ~A" i j k)))   

;;do*: a utiliser pour des initialisations imbriquées

;;*2: BLOCK
;;(un progn dont on peut sortir brutalement)
(block toto
  (setf x 1)
  (incf x)
  (return-from toto x)
  (format t "never see this"))
;;2

;;*3: DOLIST

(setf l '())
(dolist (x '(1 2 3) (nreverse l))
  (format t " ~A:" x)
  (push x l))
;;1: 2: 3:
;;(1 2 3)

;;*4: DOTIMES
;;(n iterations de 0 à n-1)

(setf l '())
(dotimes (i 5 l)
  (push i l))
;;(4 3 2 1 0)

;;*5: LOOP

(defun iota (n)
  (loop for i from 1 to n collect i))  

(iota 5)
;;(1 2 3 4 5)

(loop for n from 1 to 10
      when (oddp n)
      collect n)
;;(1 3 5 7 9)


(loop for i from 0
      for j downfrom 5
      for k = (+ i 1)
      while (< i 5)
      collect k)
;;(1 2 3 4 5)   

(loop for i from 0
      while (< i 5)
      do (format t "ff"))
;;ffffffffff

(loop for x = 5 
      for y from 0 to 5
      collect (+ x y))
;;(5 6 7 8 9)

;;cf ACL p242-3 for other loop examples


;;;-----------------------------------------------------------------------------
;;; *Les iterateurs MAP/APPLY
;;;-----------------------------------------------------------------------------

;;*1: mapcar & co
(mapcar #'1+ '(1 2 3))
;;(2 3 4)


(mapcan #'list
	'(a b c)
	'(d e f))
;;(A D B E C F)
;;au lieu de ((A D) (B E) (C F))

;;splicing
(mapcan #'identity '((1 2 (3 4)) (5 6 7)))
(1 2 (3 4) 5 6 7)


;;mapc ne crée pas de liste résultat
;;il retourne la première liste passée en argument
;;et se contente de parcourir

(mapc #'list
      '(a b c)
      '(d e f))
;;(A B C)

(mapc #'(lambda (x y) (format t "~A ~A / " x y))
      '(a b c) '(d e f))
;;* A D / B E / C F / 
;;(A B C)

;;pour les sequences: (map nil etc... )


;;*2: map-into manipule listes/tableaux/...

(setf a (list 1 2 3 4))
(setf b (make-array 4 :initial-contents '(10 10 10 10)))

;;on affecte la variable a à l'addition des vecteurs a et b
(map-into a #'+ a b)
;;(11 12 13 14) [c'est a]

;;equivalent:
(map-into a #'(lambda (x y) (+ x y)) a b)

(map-into b #'(lambda (x) (+ 1 x))
	  a)
;;#(2 3 4 5) [c'est b, l'array]


;;un tableau resultat existant  etant necessaire, voici une methode pour ne pas
;;modifier la source etne pas allouer de tableau
(map-into a #'(lambda (x) (progn (format t "~A " x) x))
	  a)


(apply #'+ '(1 2 3))
;;6


(map 'list #'- '(1 2 3 4))
;;(-1 -2 -3 -4)
;;#\: specificateur de char 
(map 'string #'(lambda (x) (if (oddp x) #\1 #\0))
     '(1 2 3 4))
;;"1010"


(map 'list #'(lambda (x) (format t "~A " x)) a)
;;11 11 11 11 :affichage du vecteur a
;;(NIL NIL NIL NIL) :retourne (nil: par format)


;;**pour ne pas etre embete avec les types de sequecncs!
(map 'list #'1+ (make-array 5 :initial-contents (iota 5)))
;;(2 3 4 5 6)


;;;-----------------------------------------------------------------------------
;;; *Conditions
;;;-----------------------------------------------------------------------------

(setf x 5)

;;COND
(cond ((= a 1) (setq a 2))
      ((= a 2) (setq a 3))
      ((and (= a 3) (floor a 2)))
      (t (floor a 3)))


;;IF prend 1 condition et 2 expressions (si cond=t) et (sinon)
(if (> x 2)
    (format t "x>2")
    (format t "x<=2"))
;;x>2
;;NIL car on renvoie rien avec les format


;;WHEN execute *une serie* d'instructions SSI la conditions est vraie
(when (> x 2)
  (format t "x>2")
  (format t "~%x n'est pas <=2")
  t)
;;x>2
;;x n'est pas <=2
;;T
;;resultat avec (setf x 1)
;;NIL


;;UNLESS execute la série d'instruction SAUF SI la condition est vraie

(unless (< x 0)
  (format t "x>=0")
  t)
;;x>=0
;;T



;;CASE - i.e switch 
;;si la clause otherwise n'est pas presente, si aucun pattern n'est matché
;;alors on renvoie NIL 
(dolist (k '(0 1 2 3 :four #\v () t 'other "vectra" nil))
  (format t "~S "
	  (case k
	    ((0 1 2) "0, 1 ou 2")
	    (3 "3")
	    (nil "pas de (val a comparer), jamais matché")
	    ((nil) "nil")
	    ((:four #\v) "4/v")
	    ((t) "t")
	    (("vectra") "project")
	    (otherwise "tout le reste"))))

;;reponses:
;;"0, 1 ou 2" "0, 1 ou 2" "0, 1 ou 2" "3" "4/v" "4/v"
;;"nil" "t" "tout le reste" "tout le reste" "nil" 
;;*****Rem: marche pas avec les strings*****

;;ECASE, CCASE: pas d'otherwise, renvoient une erreur si pas de match 


;;;-----------------------------------------------------------------------------
;;; *valeurs multiples
;;;-----------------------------------------------------------------------------


(values 1 2 3)
;;1
;;2
;;3
;;* 

;;on ne recupere pas comme ca des valeurs multiples
(setf x (values 1 2 3))
x
;;1

;;la bonne methode
(multiple-value-bind (x y z) (values 1 2 3)
  (list x y z 4))
;;(1 2 3 4)

;;si il n'y a pas assez de valeurs, NIL
(multiple-value-bind (x y z) (values 1 2 )
  (list x y z 4))
;;(1 2 NIL 4)

(multiple-value-call #'+ (values 1 2 3))
;;6

(multiple-value-list (values 'a 'b 'c))
;;(A B C)

;;intérêt: pouvoir traiter les valeurs de retour d'une fction si leur nombre
;;varie ou est inconnu
;;(pour mvalue-bind, il faut connaitre le nombre a l'avance)


;;;-----------------------------------------------------------------------------
;;; *Fermetures 
;;;-----------------------------------------------------------------------------


;;exemple d'utilisation de fermetures
(defun creer-compteur ()
  (let ((cpt -1))
    (lambda ()
      (progn (incf cpt)
	     cpt))))

(setf compteur (creer-compteur))
(funcall compteur)
;;0
(funcall compteur)
;;1

;;expl: creer-compteur renvoie une lambda + son environnement de creation qui est conserve.
;;l'envt en question etait celui d'un let, donc local. un appel de compteur ne fait que modifier cet envt.


;;;-----------------------------------------------------------------------------
;;; *Fichiers
;;;-----------------------------------------------------------------------------

;;keywds file 

;; 1-version directe

(setf path (make-pathname :name "./box.asc"))
(setf stream (open path :direction :input))
(read-line stream)
-->"Ambient light color: Red=0.039216 Green=0.039216 Blue=0.039216"
(read-line stream)
-->""
(read-line stream)
-->"Named object: \"Object01\""
(close stream)

;;une erreur survient si on atteint la fin du fichier

;; 2- assure que le fichier sera ferme 
;;meme si le corps plante ou quitte
(with-open-file (stream path :direction :input)
  (dotimes (i 4)
    (read-line stream)))

;;plus besoin de fermer le fichier


;;tests fin de fichier
(setf path (make-pathname :name "./box.asc"))
(setf s (open path :direction :input))
(do ((l (read-line s) (read-line s nil 'eof)))
    ((eq l 'eof) 'GLOS-END-OF-FILE-REACHED)
  (format t "~%~A" l))
(close s)

(let ((path (make-pathname :name "toto.txt")))
  (with-open-file (s path :direction :output :if-exists :supersede)
    (dotimes (i 5 'nil)
      (write-line (mkstr "::" i) s))))
;;::0
;;::1
;;::2
;;::3
;;::4


;; recuperer une liste de (paths de) fichiers correspondant a un motif
(setf l (directory "*.txt"))
;;(#p"/users/these/bailleul/These/Nurbs/HippocampusR.txt"
;; #p"/users/these/bailleul/These/Nurbs/LCaudateTset.txt" ...)

(pathname-name (car l))
;;HippocampusR
(pathname-type (car l))
;;txt
;;;; for other functions: http://dept-info.labri.u-bordeaux.fr/~strandh/Teaching/Programmation-Symbolique/Common/HyperSpec/Body/sec_19-1-2.html


(with-open-file (stream (make-pathname :name input-file-name)
			:direction :input)
  (loop for line = (read-line stream nil) ;; renvoie nil en fin de stream
	while line                        ;; tt que != nil
	collect (str-tokenize line)))


;;;-----------------------------------------------------------------------------
;;; Les erreurs
;;;-----------------------------------------------------------------------------

;;* 1: Unwind-Protect

;;prend en argument 1 expression speciale, puis n expressions
;;evalue la premiere, puis toutes les autres meme si la premiere plante
;;retourne la valeur de la 1ere exp.

(catch 'abort 
  (unwind-protect 
      (throw 'abort 99)
    (format t "a")
    (format t "b")))
;;ab
;;99

;;attention: si dans les n expressions restantes il y en a une qui plante,
;;leur evaluation est stoppee
(setf x 1)
(catch 'abort 
  (unwind-protect 
      (throw 'abort 98)
    (throw 'abort 99)
    (setf x 5)))
;;99
x
;;1

;;recapitulation: on n'annulle pas l'erreur, on la laisse survenir mais on
;;fait des traitements (corrects) avant.
(unwind-protect
    (push 1 nil)
  (setf x 1)
  (setf y 2))
;;erreur, on peut pas pusher 1 dans nil
x
;;1
y
;;2


;;* 2: ignore-errors 

(setf x 0 y 0)   
(ignore-errors
  (push 1 nil)
  (setf x 1)
  (setf y 2))
;;l'erreur s'affiche, mais n'interrompt pas le debugger (x et y =0 tjrs)

;;exemple de l'hyperspec
(defun load-init-file (program)
  (let ((win nil))
    (ignore-errors			;if this fails, don't enter debugger
      (load (merge-pathnames (make-pathname :name program :type :lisp)
			     (user-homedir-pathname)))
      (setq win t))
    (unless win (format t "~&Init file failed to load.~%"))
    win))
 
(load-init-file "no-such-program")
;;>> Init file failed to load.
;;NIL



;;;-----------------------------------------------------------------------------
;;; *Strings
;;;-----------------------------------------------------------------------------


(setf s (make-array 8
		    :element-type 'character
		    ;;		    :fill-pointer 3 
		    :initial-element #\ 
		    :adjustable t ))
;;"aaa"

(char s 2)
;;#\a


(setf s " \"rre\"rr tt ")
(str-next-token s) ;IOstream.lisp

(setf s2 "eee fffd ff\" X:8.5 ddd")
(str-fetch-arg "X:" s2) ;IOstream.lisp


;;;*with-output-to-string

;;creer une 'growable' string
(setq fstr (make-array '(0) :element-type 'base-char
		       :fill-pointer 0 :adjustable t)) 

(with-output-to-string (s fstr)
		       (format s " ~A" 1)
		       (format s " ~A" 2))    
;;fstr
;;" 1 2"
;; le vector se comporte comme une file!!!


;; Position d'un element (sequence)
(position (vec3-make 4.0 4.0 4.0)
	  (vec3-parse-line "5 6 9 7 8 9 4 4 4 5")
	  :test #'vec3-equal)
;;2


(concatenate 'string "aa " "bb" "cc")
;;"aa bbcc"

(mkstr "aa is:" 1 "and bb is:" 2)
;;"aa is:1and bb is:2"


;;;-----------------------------------------------------------------------------
;;; Arrays
;;;-----------------------------------------------------------------------------


(setf g (make-array 5 :initial-element 0 :adjustable t :fill-pointer t))
(setf h g)
(setf (aref g 0) 9.99)
g
#(9.99 0 0 0 0)
h
#(9.99 0 0 0 0)




;; Un tableau comme une file redimensionnable
(setf v (make-array 10 :fill-pointer 0
		    :initial-element nil :adjustable t))

v
;;#()
(dotimes (i 10 v)
  (vector-push i v))
;;#(0 1 2 3 4 5 6 7 8 9)
(vector-push 10 v)
;;NIL (vecteur trop petit)

(adjust-array v (+ (length v) 10) :initial-element nil)
(dotimes (i 9 v)
  (vector-push i v))
;;#(0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8)

(vector-push 9 v)
;;19 position du dernier element pushe
(vector-push 9 v)
;;NIL (vecteur trop petit)

(vector-pop v)
;;9
(vector-pop v)
(vector-pop v)
;;8 7
v
;;#(0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6)
;;Les elements poppes sont encore la!
(aref v 19)
;;9
(vector-push 'niou v)
;;17
v
;;#(0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 NIOU)
(fill-pointer v) position du curseur d'insertion par vector-push
;;18

;;Accesseur generique
(aref v 2)
;;1

;;Accesseur optimise pour simple-vector
(svref v 2) ;marche pas pour v!
(type-of v)
;;(VECTOR T 20)
(type-of (make-array 5 :initial-element (iota 5)))
;;(SIMPLE-VECTOR 5)



;; push "securise" d'elements: le vecteur est etendu d'un nombre
;; d'elements optionnel (x? par defaut), et fixes a 0 sauf le premier bien sur.
(vector-push-extend 1 v)  ;; on etend de x elements. ils sont fixes a 0,
;;on peut y acceder par aref, mais ils sont pas affiches et pas pris en ligne de compte par length.
;;on peut faire (setf (fill-pointer v) (+ (length v) 3)) pour y "remedier"
(vector-push-extend 1 v 1) ;;on etend d'un seul element! 


;;;-----------------------------------------------------------------------------
;;; Types et conversions
;;;-----------------------------------------------------------------------------

(type-of 2.3)
;;SINGLE-FLOAT
(typep 2.3 'single-float)
;;T

(coerce 2 'single-float)
;;2.0
(coerce "2.3" 'single-float)
;;erreur: peut pas avec string


;;;-----------------------------------------------------------------------------
;;; Curryfication
;;;-----------------------------------------------------------------------------

;;ACL
(defun curry (fn &rest args)
  #'(lambda (&rest args2)
      (apply fn (append args args2))))

(funcall (curry #'+ 3) 1 1)
;;5

(defun myplus (x y z)
  (+ x y z))
(funcall (curry #'myplus 1 1) 1)
;;3

;;usage de la macro lcurry (glos)
(setf myfun (lambda (x y z)
	      (format t "~%x:~A y:~A z:~A" x y z)))
(funcall (lcurry myfun 1 nil 3) 2)
;; x:1 y:2 z:3

;;recuperer du lambda-code a partir d'une fction
(function-lambda-expression #'myinc)
;;(LAMBDA (X) (BLOCK MYINC (+ 1 X)))
;;NIL
;;MYINC

(defun add2 (x y)
  (+ x y))

(funcall (lcurry add2 (car x) (cdr x)) '(1 2))

(add2 (list-values '(1 2)))


;;;-----------------------------------------------------------------------------
;;; *Listes
;;;-----------------------------------------------------------------------------


(setq test1 (list 1 3 4 6 7))
(setq test2 (list 2 5 8))
(merge 'list test1 test2 #'<)
;;=>  (1 2 3 4 5 6 7 8)


;;2 formes equivalentes 
(setf l (list 1 2 3 4 5))
(setf l '(1 2 3 4 5))
;;(1 2 3 4 5)

;;quote, backquote, splicing

(setf l '(1 2 3 4 5))
;;butlast: pas d'effet de bord sur l
(butlast l)
;;(1 2 3 4)
;;l = (1 2 3 4 5) encore

(some #'identity '(nil nil nil 1 nil))
;;1


(setf l (iota 5))
(append l (list 6 7))
;;(1 2 3 4 5 6 7) l inchangee: append est fonctionnel
(nconc l (list 6 7))
;; l changee!!
;;** attention: pb si l est = a (list)**

(setf l (iota 5))
(setf tete (pop l))


;;;; Listes d'association.
(setf alist (list (cons 1 "un") (cons 2 "deux") (cons 3 "Troie")))
;;((1 . "un") (2 . "deux") (3 . "Troie"))
(assoc 1 alist)
;;(1 . "un")

(member v2 vl :test #'(lambda (x y) (vec3-equal x y)))


;;select sublists
(setf l1 '(0 1 2 3 4 5 6))
(subseq l1 3)
;;(3 4 5 6)
(subseq l1 3 5)
;;(3 4)


;;fonctionnel
(remove 0 '(0 1 2 0))
(remove-if #'(lambda (x) (= 0 x))  '(0 1 2 0))
;;(1 2)

;;also:
(append '(1 2 3) '(4 5)) ;;(1 2 3 4 5)


;; count occurences (hyperspec)
(count 1 '(1 2 1 1)) ;;3
(count #\a "how many A's are there in here?");; =>  2 a
(count-if-not #'oddp '((1) (2) (3) (4)) :key #'car) ;;=>  2 evens: 2 & 4 
(count-if #'upper-case-p "The Crying of Lot 49" :start 4) ;; =>  2 (C & L)


;;cf also sort on this page

;;;-----------------------------------------------------------------------------
;;; Profiling
;;;-----------------------------------------------------------------------------


(load "library:subsystems/defsystem-library.x86f")
(use-package :profile)

;;activer on/off le profiling d'une fonction
(profile asc-read-object)
(unprofile asc-read-object)

(asc-parse-file "../asc/dolphins2.asc")

;;consulter les ressources consommees par une fction
(report-time asc-read-object)
;;mettre a zero ce total (cumulatif)
(reset-time asc-read-object)



;;;-----------------------------------------------------------------------------
;;; *Flet
;;;-----------------------------------------------------------------------------

(flet ((myinc (x) (+ 1 x)))
  (myinc 5))
;;6
;;rem: la portee de myinc n'excede pas le flet 


(flet ((inc1 (x) (+ 1 x))
       (inc2 (y) (+ 1 (inc1 y))))
  (inc1 5)
  (inc2 6))
;;erreur, on peut pas utiliser une definition precedente

(labels ((inc1 (x) (+ 1 x))
	 (inc2 (y) (+ 1 (inc1 y))))
  (inc1 5)
  (inc2 6))
;;8
;;avec labels, on peut! On peut meme definir des fonctions recursives.


(prog (var1 var2 (var3 init-form-3) var4 (var5 init-form-5))
      declaration*
      statement1
      tag1
      statement2
      statement3
      statement4
      tag2
      statement5
      ...
      )


;; labels/flet peut servir creer un environnement ou on pourra definir
;; 1 ou plusieurs fonctions qui garderont les fonctions auxilliaires
;; dans leur closure.
(labels ((inc (x) (1+ x)))
  (defun toto (x)
    (+ 2 (inc x))))


;;;-----------------------------------------------------------------------------
;;; *CLOS
;;;-----------------------------------------------------------------------------


(defclass point ()
  (x y z))

(setf p1 (make-instance 'point))

(setf (slot-value p1 'x) 5.0)
(slot-value p1 'x)
;;5.0
;;erreur si on accede a un slot non-affecté

;;!!! cette definition ne supporte pas de commentaire ""
(defclass point ()
  ((x :accessor access-x);;r/w
   (y :reader get-y)
   (z :writer set-z)))

(setf p2 (make-instance 'point))

(setf (access-x p2) 5.0)
(access-x p2)
;;5.0

;;attetion a l'ordre des arguments
(set-z 8.3 p2)


(defclass 2cv ()
  ((couleur :accessor get-couleur
	    :initarg :couleur;;mot-clé argument d'instanciation
	    :initform 'bleue);;valeur par défaut 
   (etat :accessor get-etat
	 :initform '(délabré et sale))
   (surnom :initform '22che 
	   :allocation :class)))

;; :allocation :class = variable de classe
;; :initarg :couleur = mot-clé argument d'instanciation
;; :initform 'bleue = valeur par défaut

(setf vectra (make-instance '2cv :couleur 'bleucrade)) 
(get-etat vectra)
;; (|DéLABRé| ET SALE)
(get-couleur vectra)
;; BLEUCRADE

(slot-value vectra 'surnom)
;; 22CHE

(slot-value (make-instance '2cv :couleur 'rouge) 'surnom)
;; 22CHE

;;remarque: si on change cette valeur sur une instance, la valeur est modifiée pour ttes les autres 

;;remarque: tu changes la definition de la classe, et les instances deja créées sont mises a jour



(defclass vehicule ()
  ((couleur :accessor get-couleur
	    :initarg :couleur
	    :initform 'rouge)
   (etat :accessor get-etat
	 :initform '(normal))))


(defclass objet-en-metal ()) 

;;la classe voiture hérite (attributs, méthodes) des classes vehicule/objet-en-metal
(defclass voiture (vehicule objet-en-metal)
  ((surnom :accessor access-surnom
	   :initform 'titine 
	   :allocation :class)))

;;ajout d'une méthode à la classe voiture
(defmethod demarrer ((v voiture) commentaire)
  (format t "~% ~A demarre et son conducteur dit: ~A"
	  (access-surnom v)
	  commentaire)
  t)


(setf 2cv (make-instance 'voiture :couleur 'bleu-crade))

(demarrer 2cv "allez demarre")
;;* TITINE demarre et son conducteur dit: allez demarre
;;T

;;Ajouter des actions à effectuer avant ou après l'exécution d'une méthode

(defmethod demarrer :before (v voiture)
  (let ((surnom (access-surnom v)))
    (format t "~% ~A: tchitchitchi" surnom)))

(defmethod demarrer :after (v voiture)
  (let ((surnom (access-surnom v)))
    (format t "~% ~A: Vrouuuuumm reup reup reup!" surnom)))

(demarrer 2cv "allez demarre")
;;TITINE: tchitchitchi
;;TITINE demarre et son conducteur dit: allez demarre
;;TITINE: Vrouuuuumm reup reup reup!


;;methodes virtuelles:

;;pas de code à la déclaration de la méthode virtuelle
(defgeneric avance ((v vehicule))
  (:documentation "avance le vehicule"))

(defmethod avance ((v vehicule))
  :documentation "méthode "avance" par défaut pour tout véhicule"
  (format t "~% vroum le véhicule"))
  
(avance 2cv)
;;* vroum le véhicule
;;NIL

(defmethod avance ((v voiture))
  :documentation "méthode pour les voitures seulement"
  (format t "~% vroum la tuture"))

(avance 2cv)
;;* vroum la tuture
;;NIL


;;initialization
;;attention, si on lit des propriétés non affectées, erreur

(defclass trucs ()
  ((val1 :initform 1
	 :accessor acc1)
   (val2 :accessor acc2)))

(defmethod initialize-instance :after ((tr trucs) &key)
  (format t "~% val1: ~A" (acc1 tr))
  (setf (acc1 tr) 5)
  (setf (acc2 tr) 6)
  (format t "~% val1: ~A ~% val2: ~A" (acc1 tr) (acc2 tr)))

(setf t1 (make-instance 'trucs))
;;* val1: 1
;; val1: 5 
;; val2: 6
;;#<TRUCS {48050C1D}>


;;;-----------------------------------------------------------------------------
;;; *Sauts non locaux
;;;-----------------------------------------------------------------------------


(defun f (x)
  (throw 'hello 345))

(defun g (a)
  (catch 'hello
    (format t "~% avant")
    (f a)
    (format t "~% apres")))

(g 356)
;; avant
;;345


;;;-----------------------------------------------------------------------------
;;; *Maths
;;;-----------------------------------------------------------------------------


(incf x)
(decf x)
;;(inc/dec)rémentent des variables, et non des valeurs

;;1+/1- plus généraux
(1+ x)
(1- x)
(1+ 5)


(ceiling 2.4)
;; 3
;;-0.5999999 : 2nd return value

;;cf: floor, round

(sort '(1 6 5 5 4 7 1 2) #'<)
;;(1 1 2 4 5 5 6 7)


;;;-----------------------------------------------------------------------------
;;; *Compilation - debogueur
;;;-----------------------------------------------------------------------------
;; present dans ~/.cmucl-init

(defun declaim-fast ()
  "set compiler options to fast not-debuggable code"
  (declaim (optimize (speed 3)
		     (compilation-speed 0)
		     (safety 0)
		     (debug 0)))) 


(defun declaim-easy ()
  "set compiler options to slow reliable debuggable code"
  (declaim (optimize (speed 1)
		     (compilation-speed 0)
		     (safety 2)
		     (debug 3))))


(defun declaim-debug ()
  "set compiler options to slow reliable debuggable code"
  (declaim (optimize (speed 0)
		     (compilation-speed 0)
		     (safety 3)
		     (debug 3))))



(compile-file "./myfile.lisp") ;; marche aussi avec un path
(load  "./myfile.sparcf") ;ou .x86f sur linux:
;;necessaire pour que les fonctions compilees soient employees
;;#'findspan -> #<Interpreted Function FINDSPAN {4004FA49}> si non compilee
;;#'findspan -> #<Function FINDSPAN {40B9CBC1}> si compilee

(compile 'fun_name)



;;declarer la fonction
(defun funtest (x) (+ 1 x))
(compile 'funtest)
(break)
;; on passe en mode debugeur
;;0] ll #'funtest
;;0: #'(LAMBDA (X) (BLOCK FUNTEST (+ 1 X)))
;;1: (+ 1 X)
;;0] br 1
;;(+ 1 X)
;;1: 1 in FUNTEST
;;Added.
;;0]
;;retour au top-level et exec de la fonction pour activer le breakpt

;;step n ou seul
;;delete-breakpoint n ou seul pour tous

;;Tres souvent: descendre a repetition
;;down 
;;pour voir le code
;;source
;;source nligne
;;pour voir les variables
;;list-locals

;; TODO : qd un programme plante:
;; source
;; changer de frame (utile pour savoir ou un programme a coince)
;;f 1  -- ou down

;;interroger cmucl sur les plus grosses variables presentes:
(room t)

;;afficher arguments d'appels de chaque fonction
(trace function-name)

;;Tres utile: inserer (break) dans le code de la fonction: passage en mode interactif
;;on peut afficher la valeur des variables et faire des (setf keeped variable)


;;Ne pas oublier de baisser le niveau d'optimisation pour debugguer!
;;evt, loader le fichier source, et verifier par #'fun que notre
;;fonction est bien en mode interperte.

;;** gros avtage par rapport au debugguer C: on peut faire des
;;affectations des valeurs locales trouvees qui resteront apres debog
;;on peut lancer des procedures a partir du degogueur sans coredump
;;quand on passe a l'objet, c'est toujours lisible


;;Desassemblage: voir le keude ASM
(defun int-add1 (n)
  (declare (fixnum n)
	   (optimize (speed 3) (safety 0) (debug 0)))
  (the fixnum (1+ n)))
(disassemble 'int-add1)


;;;-----------------------------------------------------------------------------
;;; *Macros
;;;-----------------------------------------------------------------------------


(defmacro ntimes (n &rest body)
  `(do ((x 0 (+ x 1)))
       ((>= x n) ,n)
     ,@body))
    
(symb 'toto_ 123) ;; => TOTO_123 symbole!


;;;-----------------------------------------------------------------------------
;;; tables de hachage 
;;;-----------------------------------------------------------------------------


(setf ht (make-hash-table))
;;#<EQL hash table, 0 entries {803EAB6D}>

(setf (gethash 'key1 ht) "key1 value")
(gethash 'key1 ht)
;;"key1 value"
;;T

;;*** attention ***: eql est employee pour comparer les cles: marche
;;avec les symboles, keywords, constantes, mais pas les chaines.
;;idee: (symb "machainecle")

(setq table (make-hash-table :test #'equal))
;;#<HASH-TABLE EQUAL 0/139 46145547>

(setf (gethash "one" table) 1)
;;1
(gethash "one" table)
;;1
;;T

;; iteration:

(setq t2 (make-hash-table :test 'equal))

(dotimes (i 10)
  (setf (gethash i t2) i)) 

(maphash #'(lambda (key val)
	     (format t "~%~A ~A" key val))
	 t2)
;;0 0
;;1 1
;;...
;;9 9
;;NIL
;;rappel: on peut remplacer la lambda par un nom de fonction

;;cf aussi: WITH-HASH-TABLE-ITERATOR 

;; vider une table:
(clrhash t2)


(flet ((ite (x y) (format t "~%~A ~A" x y)))
  (maphash #'ite t2))



;;;-----------------------------------------------------------------------------
;;; Types
;;;-----------------------------------------------------------------------------


(typep  13134654654 'integer)
;;T
(type-of 1313465465)
;;BIGNUM


(coerce (read-from-string str)
	  'real)


;;;-----------------------------------------------------------------------------
;;; Garbage collection - maybe CMUCL specific
;;;-----------------------------------------------------------------------------
;; http://www.cons.org/cmucl/doc/gc-tuning.html

(room)
;;gc memory usage


(format t "~D" ext:*bytes-consed-between-gcs*)

(let ((ext:*bytes-consed-between-gcs*
       (* ext:*bytes-consed-between-gcs* 5))
      (acc (list)))
  (format t "~D" ext:*bytes-consed-between-gcs*)
  (dotimes (i 1000 nil)
    (setf acc (cons 1 acc))))

(format t "~D" ext:*bytes-consed-between-gcs*)

;;8000000
;;NIL
;;* 40000000
;;NIL
;;* 8000000


;;;-----------------------------------------------------------------------------
;;; Constants, global variables
;;;-----------------------------------------------------------------------------
;;D Lamkins

;;You should define global constants using DEFCONSTANT. From the
;;viewpoint of reading a Lisp program, the distinction between
;;DEFPARAMETER and DEFCONSTANT is that the value defined by
;;DEFPARAMETER could concievably be altered by the user after the
;;program is compiled, but a DEFCONSTANT value will never change. A
;;good Lisp compiler will take advantage of DEFCONSTANT declarations
;;to perform classical optimizations such as constant folding or
;;compiling immediate load instructions.

;;Fewer Lisp programmers follow a naming convention for constants. The one I use puts a leading and trailing plus sign on the name of the constant, as in +RTS-OPCODE+. 

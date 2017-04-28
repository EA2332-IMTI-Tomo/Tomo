;;; ----------------------------------------------------------------------------
;;; $RCSfile: Tools.lisp,v $
;;; $Author: Vectra Project $
;;; $Date: 2000/10/23 06:13:39 $
;;; $Revision: 1.1.1.1 $
;;; Extended version taken from GLOS project
;;; ----------------------------------------------------------------------------


;;(declaim (optimize (speed 3)
;;		   (compilation-speed 0)
;;		   (safety 0)
;;		   (debug 0)))


;;;-----------------------------------------------------------------------------
;;; Useful Macros
;;;-----------------------------------------------------------------------------


;; sets given global variable to new value inside body
(defmacro with-global-var-excursion (var-name var-tempval &body body)
  (let ((var-previousval (gensym)))
    `(let ((,var-previousval ,var-name)) 
       (progn (setf ,var-name ,var-tempval)
	      ,@body
	      (setf ,var-name ,var-previousval))))) 


;;Paul Graham's 'Ansi Common Lisp'
(defmacro while (test &body body)
  `(do ()
       ((not ,test))
     ,@body))


(defmacro dorange (iterator start end rvalue &body body)
  "executes body for iterator values going from start to end both
included, in/decreasing depending on start>/<end."
  (let ((inc-fun (if (> end start) #'1+ #'1-))
	(comp-fun (if (> end start) #'> #'<)))
    `(do* ((,iterator ,start (funcall ,inc-fun ,iterator)))
	 ((funcall ,comp-fun ,iterator ,end) ,rvalue)
       ,@body)))


(defmacro progt (&body exps)
  "execute given expressions in sequence and returns T at the end."
  `(progn ,@exps
	  t))


;;Example: (setf myfun (lambda (x y z)
;;	      (format t "~%x:~A y:~A z:~A" x y z)))
;;(funcall (lcurry myfun 1 nil 3) 2)
;; --> x:1 y:2 z:3
(defmacro lcurry (lambdaform &body args)  
  "returns a curried version of given lambda"
  (let ((myargs '())
	(symlist '()))
    (progn (dolist (x args)
	     (if x
		 (push x myargs)
		 (let ((sym (gensym)))
		   (push sym symlist)
		   (push sym myargs))))
	   `(lambda ,(nreverse symlist)
	      (funcall ,lambdaform ,@(nreverse myargs))))))


;;consider multiple-value-setq for multiple-value assignment
(defmacro nsetf (value &body vars)
  "declare n variables initialized to an unique value."
  (let ((acc '()))
    (dolist (x vars)
      (push x acc)
      (push value acc))
    `(setf ,@(nreverse acc))))


(defmacro list-match-setf (lsrc ldest)
  "setf n variables to n values, all included in separate lists"
  (let ((gencode (lambda (x y) `(setf ,x ,y))))
    `(progn ,@(mapcar gencode lsrc ldest)
	    t)))


(defmacro list-match-incf (lsrc ldest)
  "incf n variables of n increments, all included in separate lists"
  (let ((gencode (lambda (x y) `(incf ,x ,y))))
    `(progn ,@(mapcar gencode lsrc ldest)
	    t)))


;;impossible to do (mapcar #'mymacro-1arg '(argval1 argval2 ... argvaln)
(defmacro map-macro (macroname arglist)
  "calls macro on every successive value of arglist"
  (let ((acc (list)))
    (dolist (arg arglist
		 (cons 'progn (nreverse (cons t acc))))
      (push `(,macroname ,arg) acc)))) 


;;the only way to implement scheme's (define fun2 fun1)
;;cf on-lisp p 215
(defmacro bindfun (fun1 fun2)
  "binds a name function to another"
  `(defmacro ,fun1 (&rest body)
     `(,',fun2 ,@body)))


;;;-----------------------------------------------------------------------------
;;; Useful functions from Paul Graham
;;;-----------------------------------------------------------------------------



;;Paul Graham's 'Ansi Common Lisp'
(defun curry (fn &rest args)
  #'(lambda (&rest args2)
      (apply fn (append args args2))))


;;Paul Graham's 'On Lisp'
(defun mkstr (&rest args)
  (with-output-to-string (s)
			 (dolist (a args) (princ a s))))


;;Paul Graham's 'On Lisp'
(defun symb (&rest args)
  (values (intern (apply #'mkstr args))))


;;;-----------------------------------------------------------------------------
;;; GLOS Reader Macros
;;;-----------------------------------------------------------------------------

;; if 2 reader macros inverted, matrix one fails
;; reader macro for accessing vector slots

;;both cannot coexist

;;(progn
;;  (set-macro-character
;;   #\} (get-macro-character #\)))
;;  
;;  (set-dispatch-macro-character
;;   #\# #\{
;;   #'(lambda (stream char1 char2)
;;       (declare (ignore char1 char2))
;;       (let* ((trike (read-delimited-list #\} stream t))
;;	      (v (car trike))
;;	      (i (cadr trike)))
;;	 `(the single-float (aref ,v ,i))))))



;; reader macro for accessing mat44 elements
;; #{m 2 3} -> 1.0
(progn
  (set-macro-character
   #\} (get-macro-character #\)))

  (set-dispatch-macro-character
   #\# #\{
   #'(lambda (stream char1 char2)
       (declare (ignore char1 char2))
       (let* ((trike (read-delimited-list #\} stream t))
	      (m (car trike))
	      (i (cadr trike))
	      (j (caddr trike)))
	 (if (numberp i) 
	     (if (numberp j)
		 (let ((index (the fixnum (+ (coerce j 'fixnum) (* (coerce i 'fixnum) 4)))))
		   `(the single-float (aref ,m ,index)))
					; j n'est pas un nombre
		 (let ((index (the fixnum (* (coerce i 'fixnum) 4))))
		   `(the single-float (aref ,m (+ ,j ,index)))))
	     `(the single-float (aref ,m (+ ,j (the fixnum (* ,i 4))))))))))




;;;-----------------------------------------------------------------------------
;;; Useful list utilities
;;;-----------------------------------------------------------------------------


;;arguments declaration & update observe the "do" syntax
;;Example: (generate-list 5 (lambda (x y) (+ x y))
;;                        (a 0 (+ a 2)) (b 1 (- b 1)))
;; --> (1 2 3 4 5)
(defmacro generate-list (len genfunc &body paramlist)
  "generate a list of "len" elements, result of the applying the given
lambda to a set of possibly updated arguments"
  (let ((params (mapcar #'car paramlist))
	(k (gensym))
	(acc (gensym)))
    `(do* (,@paramlist
	   (,k 0 (1+ ,k))
	   (,acc '()))
	 ((>= ,k ,len) (reverse ,acc))
       (push (funcall #',genfunc ,@params) ,acc))))


(defmacro iota (n)
  "creates a list of n integers from 1 to n"
  `(generate-list ,n identity (x 1 (1+ x))))


(defun iota-range (start end)
  "creates a list of n integers from start to end included "
  (let* ((len (1+ (- end start))))
    (generate-list len identity (x start (1+ x)))))


(defun iota-from (start length)
  "create a list of length integers from start included"
  (iota-range start (+ start (1- length))))
  

(defun iota-centered (center half-length)
  (iota-range (- center half-length) (+ center half-length)))


(defmacro unilist (n elt)
  "creates a list of n elements"
  `(generate-list ,n identity (x ,elt)))


(defmacro remove! (item l)
  "removes all occurences of item in given list (side-effect)"
  `(setf ,l (remove ,item ,l)))


(defun sample-list (input-list sampling-gap)
  "return a list whose entries are collected every gap+1 input-list
entries, just after input-list car. gap must be of integer type."
  (do* ((input-size (length input-list))
	(counter 0 (1+ counter))
	(output-list))
      ((= counter input-size) (nreverse output-list))
    (progn (when (zerop (mod counter (1+ sampling-gap)))
	     (push (car input-list) output-list))
	   (pop input-list))))


(defun list-values (l)
  "return values of a list"
  (apply #'values l))


;; substn '(1 a 2 b 3 c) '(a b 4 c coucou) -> '(1 2 4 3 coucou)
(defun substn (lst body)
  "multiple substitutions in a list"
  (if (null lst)
      body
      (substn (cddr lst)
	      (subst (first lst)
		     (second lst)
		     body))))


(defun splice (list-of-lists)
  "transform a list of items list into a single items list - these
unchanged items possibly including lists or whatever"
  (mapcan #'identity list-of-lists))


(defmacro append! (core-list append-list)
  "appends 2nd list to the end of the first, even nil (contrarily to
nconc)"
  `(if (endp ,core-list)
       (setf ,core-list ,append-list)
       (nconc ,core-list ,append-list)))


(defun mappopcar-vec (veclst)
  "from a list of arrays ([a b c] [A B C]...), compute a list of
arrays ([a A...] [b B...] [c C...])"
  (let ((array-lst (list)))
    (dotimes (i (length (car veclst)) (nreverse array-lst))
      (let ((result (list)))
	(mapcar #'(lambda (vec)
		    (push (aref vec i) result))
		veclst)
	(push (make-array (length result) :initial-contents (nreverse result))
	      array-lst)))))


(defun mappopcar-lst (llst)
  "from a list of list ((a b c) (A B C)...), compute a list of
arrays ([a A...] [b B...] [c C...])"
  (let ((array-lst (list)))
    (dotimes (i (length (car llst)) (nreverse array-lst))
      (let ((result (list)))
	(mapcar #'(lambda (lst)
		    (push (nth i lst) result))
		llst)
	(push (make-array (length result) :initial-contents (nreverse result))
	      array-lst)))))


(defmacro map-butlast (fun &body lists)
  "mapcar given function to the butlast of each input list"
  `(mapcar ,fun ,@(mapcar #'(lambda (x) `(butlast ,x))
			  lists)))


(defun lists-rotate (llst)
  "from a list of list ((a b c) (A B C)...), compute a list of
lists ((a A...) (b B...) (c C...))"
  (let ((list-lst (list)))
    (dotimes (i (apply #'min (mapcar #'length llst))
		(nreverse list-lst))
      (let ((result (list)))
	(mapcar #'(lambda (vec)
		    (push (nth i vec) result))
		llst)
	(push (nreverse result)
	      list-lst)))))


(defun divide-every (inlist subsize)
  "returns a list of sublists of subsize elements (from input
list). last list might contain less tahn subsize elts"
  (let ((sublists '())
	(sublists-num (ceiling (/ (length inlist) subsize))))
    (dotimes (l sublists-num (nreverse sublists))
      (push (subseq inlist (* l subsize) (* (1+ l) subsize))
	    sublists))))
      

(defun popcars (l n)
  "removes! n first list elements"
  (dotimes (i n l)
    (pop l)))


(defun poptail (l n)
  "removes! n last list elements"
  (nreverse (popcars (nreverse l) n)))


(defun shrink (l size)
  "limitates! list-size to given size, removes last element"
  (let ((tokill (- (length l) size)))
    (if (> tokill 0)
	(poptail l tokill)
	l)))


(defun firsts (n l)
  (subseq l 0 n))


(defun firsts-that (pred l)
  "return the first consecutive elements of the list that satisfy
given predicate"
  (cond ((endp l) l)
	((not (funcall pred (car l))) nil)
	(t (cons (car l)	
		 (firsts-that pred (cdr l))))))
	 

(defun lasts (n l)
  (let ((lsize (length l)))
    (subseq l (- lsize n) lsize)))


(defun pushend (elt l)
  "pushes! elt at end of list"
  (nreverse (cons elt (nreverse l))))


(defun count-from-start (val listvals)
  "returns number of occurences of value in list since its car to the
first non-equal"
  (let ((count 0))
    (dolist (x listvals count)
      (if (equal x val)
	  (incf count)
	  (return count)))))


(defun select-values (listvals listbools)
  "selects values from listvals according to successive values of
corresponding booleans list"
  (let ((acc (list)))
    (dolist (v listvals (nreverse acc))
      (let* ((arglen (length listbools))
	     (cur-arg (if (> arglen 0)
			  (pop listbools)
			  nil)))
	(when cur-arg
	  (push v acc))))))


;;(count-occurence-series 1 '(0 1 1 1 0 1 1 0 0 1 0))  => (3 2 1)
(defun count-occurence-series (val listvals)
  "return the list of consecutive counts of val in given list"
  (let ((acc (list))
	(count 0))
    (dolist (v listvals (nreverse acc))
      (if (equal v val)
	  (incf count)
	  (if (not (zerop count))
	      (progn (push count acc)
		     (setf count 0)))))))


;;position only returns pos. of 1st occurence
(defun positions (element l)
  "return a list of positions of occurences of element in given list"
  (let ((acc (list))
	(cursor 0))
    (dolist (elt l (nreverse acc))
      (when (equal elt element)
	(push cursor acc))
      (incf cursor))))
	  
		    
;;;-----------------------------------------------------------------------------
;;; Math utilities
;;;-----------------------------------------------------------------------------



        
;;;-----------------------------------------------------------------------------
;;; Array/vector iterator
;;;-----------------------------------------------------------------------------


(defun make-array-iterator (myarray &key (secure nil))
  "returns an object giving successive values of initially given
array"
  (let ((cursor -1)
	(size (length myarray)))
    (if secure
	(lambda ()
	  (progn (incf cursor)
		 (if (>= cursor size)
		     nil
		     (aref myarray cursor))))
	(lambda ()
	  (progn (incf cursor)
		 (aref myarray cursor))))))



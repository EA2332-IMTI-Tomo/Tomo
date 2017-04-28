(defun make-counter ()
  "create an auto-incremential counter"
  (let ((cpt -1))
    (lambda ()
      (progn (incf cpt)
	     cpt))))


(defmacro npush (l &body vars)
  "pushes elements to the given list"
  `(setf ,l (append ',vars ,l)))

;;(macroexpand `(npush M 7 8))
;;(SETQ M (APPEND '(7 8) M))


(defun mat44-display (m)
  (declare (type mat44 m))
  (dotimes (row 4)
    (declare (type fixnum row))
    (dotimes (col 4)
      (declare (type fixnum col))
      (format t "~F " #[m row col]))
    (format t "~%")))


(defun truep (exp)
  "true predicate"
  (not (eq nil exp)))


(defun list-or (l)
  "test if one element is not nil"
  (dolist (x l)
    (if x (return t))))


(defmacro cartesian (&body args)
  "provisoire et bâclé"
  (assert (equal args '(a b c)))
  ''((a b) (a c) (b c)))

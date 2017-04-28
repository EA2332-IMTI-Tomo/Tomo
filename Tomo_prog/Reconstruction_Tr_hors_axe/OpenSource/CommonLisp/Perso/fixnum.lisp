;; Raymond Toy
(in-package :vm)
(eval-when (:compile-toplevel :load-toplevel :execute)

(defknown fixnum+ (fixnum fixnum)
   fixnum
   (movable foldable flushable))

(define-vop (fixnum-add fast-+/fixnum=>fixnum)
   (:translate fixnum+))
)

(defun fixnum+ (x y)
   (declare (fixnum x y))
   (fixnum+ x y))

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





;;ivan boldyrev:


(EVAL-WHEN (:COMPILE-TOPLEVEL :LOAD-TOPLEVEL :EXECUTE)
	   (vm::defknown fixnum+ (fixnum fixnum)
			 fixnum
			 (vm::movable vm::foldable vm::flushable))

	   (vm::define-vop (fixnum-add vm::fast-+/fixnum=>fixnum)
			   (:translate fixnum+))

	   (vm::defknown fixnum* (fixnum fixnum)
			 fixnum
			 (vm::movable vm::foldable vm::flushable))

	   (vm::define-vop (fixnum-mul vm::fast-*/fixnum=>fixnum)
			   (:translate fixnum*))
	   )

(defun fixnum+ (x y)
  (declare (fixnum x y))
  (fixnum+ x y))

(defun fixnum* (x y)
  (declare (type fixnum x y))
  (fixnum* x y))





;;(in-package :vm)
;;(defknown fixnum+ (fixnum fixnum)
;;  fixnum
;;;;  (movable foldable flushable))

;;(define-vop (fixnum-add fast-+/fixnum=>fixnum)
;;  (:translate fixnum+))

;;(defun fixnum+ (x y)
;;  (declare (fixnum x y))
;;  (fixnum+ x y))

;;Then

;;(vm:fixnum+ most-positive-fixnum 1) => -536870912


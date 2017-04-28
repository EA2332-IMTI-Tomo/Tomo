(defun boundedp (value inf sup &key (noequal 1))
  "checks wether inf <= val <= sup. "
  (and (>= value inf)
       (<= value sup)))


(defun frame-value (value inf sup)
  "return value bounded by given borns."
  (cond ((< value inf) inf)
	((> value sup) sup)
	(t value)))


;; cf. also stats.lisp

(defun neighbors (value1 value2 precision)
  "checks wether given 2 values are close enough."
  (<= (abs (- value1 value2)) precision))


(defmacro mulf (operand mult)
  "like incf, but for multiplication"
  `(setf ,operand (* ,operand ,mult)))


(defmacro divf (operand quot)
  "like mulf, but for division"
  `(setf ,operand (/ ,operand ,quot)))


;;;;a simple syntax error like a point after string end crashes all!!
(defun function-maximum (ar1-fun &key (minimum nil)
				 (input-entries nil)
				 (start-x 0) (end-x 0) (gap 1))
  "seeks for the input value for which given function has its max(min)
value."
  (let* ((entries-list (if (eql input-entries nil)
			   (let ((entries-num (1+ (floor (/ (- end-x start-x) gap)))))
			     (generate-list entries-num identity (a start-x (+ a gap))))
			   input-entries))
	 (values-list (mapcar ar1-fun entries-list)))
    (nth (position (apply (if (eql minimum nil) #'max #'min)
			  values-list)
		   values-list)
	 entries-list)))

;;ok
;;(round (function-maximum #'(lambda (x) (* x x)) :minimum t :start-x -2 :end-x 2 :gap 0.2))
;;0


;;Maths must be loaded before compiling (image.lisp) any function using this macro
(defmacro no/0 (var expression)
  "prevents divison by zero in expression containing var as
denominator. var is returned if equals 0. Macro form allows not
evaluating expression."
  `(if (zerop ,var)
       ,var
       ,expression))


(defun sf (x)
  (coerce x 'single-float))
(defun coerce-sf (x)
  (coerce x 'single-float))


(defmacro square (exp)
  "unsecure macro sensitive to variable capture and multiple eval."
  `(* ,exp ,exp))


;;(setf a 5)
;;(square (incf a)) -> a = 7 instead 6


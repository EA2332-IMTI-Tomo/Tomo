;;; ----------------------------------------------------------------------------
;;; $RCSfile: Math3D.lisp,v $
;;; $Author: hubert $
;;; $Date: 2000/05/05 06:13:39 $
;;; $Revision: 1.1.1.1 $
;;; ----------------------------------------------------------------------------

;;GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
;;Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net)


;; ce fichier est casse pieds a compiler: si la macro
;; make-mat44-defuns et son appel ne sont pas loades avant la compil,
;; ca ne passe pas.


(defconstant +PI+ 3.141592653)
(declaim (type single-float +PI+))


(defconstant +mat44-index-list+ '((0 0) (0 1) (0 2) (0 3)
				  (1 0) (1 1) (1 2) (1 3)
				  (2 0) (2 1) (2 2) (2 3)
				  (3 0) (3 1) (3 2) (3 3)))


;;;-----------------------------------------------------------------------------
;;;-----------------------------------------------------------------------------
;; 4x4 Matrices 
;;;-----------------------------------------------------------------------------
;;;-----------------------------------------------------------------------------


(deftype mat44 () '(simple-array single-float (16)))

;; on ne peut pas utiliser le reader macro #{} ici a cause de l'ordre de la macroexpansion:
;; la reader macro est toujours expansee *avant* le reste de la macro. donc on utilise mat44-elem.


;; Constructor
(defun mat44-make-null ()
  "Create the null matrix"
  (make-array 16 :element-type 'single-float))


(defparameter *spare-matrix* (mat44-make-null))


;;;-----------------------------------------------------
;;; Matrix specific functions
;;;-----------------------------------------------------


(defmacro make-mat44-defuns ()
  "define accessors for the mat44 type"
  (let* ((index +mat44-index-list+)
	 (mat-access (lambda (mat ilist) `(#{,mat ,(car ilist) ,(cadr ilist)})))
	 (set-zero (lambda (x) `(setf ,@x 0.0)))
	 (accessors (mapcar (curry mat-access 'm) index)) ;; m existe pas!
	 (mtarget (mapcar (curry mat-access 'mtarget) index))
	 (msource (mapcar (curry mat-access 'msource) index))
	 (mat-slots (mapcar (lambda (ilist) (symb 'M (car ilist) (cadr ilist)))
			    index)))
    `(progt (defun mat44-set-null (m)
	      "Modify a matrix which become a null matrix"
	      (declare (type mat44 m))
	      (progt ,@(mapcar set-zero accessors)))
	    
	    (defun mat44-copy (mtarget msource)
	      "copy the matrix msource into mtarget"
	      (declare (type mat44 m))
	      (progt ,@(mapcar (lambda (x y) `(setf ,@x ,@y))
			       mtarget msource)))

	    (defun mat44-set (m ,@mat-slots)
	      "fill a matrix with elements"
	      (progt ,@(mapcar (lambda (x y) `(setf ,@x ,y))
			       accessors mat-slots))))))


;;effective definition of the accessors in the project
(make-mat44-defuns)


(defun mat44-identity ()
  (let ((mi (mat44-make-null)))
    (mat44-set mi
	       1.0 0.0 0.0 0.0
	       0.0 1.0 0.0 0.0
	       0.0 0.0 1.0 0.0
	       0.0 0.0 0.0 1.0)
    mi))


;;r/w accessor (aref)
;; i and j must be numbers!
(defmacro mat44-elem (m i j)
  (declare (type card8 i j))
  `(the single-float (aref ,m ,(+ (* 4 i) j))))



(defmacro row-by-col (mat-dest mat1 mat2 i j)
  (declare (type card8 i j))
  `(setf (mat44-elem ,mat-dest ,i ,j)
	 (the single-float
	   (+ (the single-float (* (mat44-elem ,mat1 ,i 0) (mat44-elem ,mat2 0 ,j)))
	      (the single-float (* (mat44-elem ,mat1 ,i 1) (mat44-elem ,mat2 1 ,j)))
	      (the single-float (* (mat44-elem ,mat1 ,i 2) (mat44-elem ,mat2 2 ,j)))
	      (the single-float (* (mat44-elem ,mat1 ,i 3) (mat44-elem ,mat2 3 ,j)))))))


(defmacro calc-line (mat-dest mat1 mat2 i)
  (declare (type card8 i))
  `(progn (row-by-col ,mat-dest ,mat1 ,mat2 ,i 0)
	  (row-by-col ,mat-dest ,mat1 ,mat2 ,i 1)
	  (row-by-col ,mat-dest ,mat1 ,mat2 ,i 2)
	  (row-by-col ,mat-dest ,mat1 ,mat2 ,i 3)))


(defconstant +mat44-identity+ (mat44-identity))
(defparameter *current-matrix* (mat44-identity))


(defun mat44-make-dup (m)
  "duplicate a matrix"
  (declare (type mat44 m))
  (copy-seq m))


(defun mat44-row (matrix rownum)
  "return a list, copy of the nth row of the given matrix"
  (generate-list 4
		 (lambda (x) #{matrix rownum x})
		 (j 0 (1+ j))))


(defun mat44-column (matrix colnum)
  "return a list, copy of the nth column of the given matrix"
  (generate-list 4
		 (lambda (x) #{matrix x colnum})
		 (i 0 (1+ i))))


(defun mat44-index-list ()
  "generate the list of indexes for the mat44 matrix type"
  (do ((i 0 (1+ i))
       (acc '()))
      ((>= i 4) (nreverse acc))
    (declare (type card8 i)) 
    (dotimes (j 4)
      (push `(,i ,j) acc))))


(defun mat44-mult (mat1 mat2 matres)
  "Multiply mat1 by mat2 and stores result into the current matrix"
  (declare (type mat44 mat1 mat2))
  (mat44-set-null matres) ;;necessary?
  (calc-line matres mat1 mat2 0)
  (calc-line matres mat1 mat2 1)
  (calc-line matres mat1 mat2 2)
  (calc-line matres mat1 mat2 3))


(defun mat44-mult-vec4 (mat vec)
  "Multiply mat by vec and raises result vec4"
  (declare (type mat44 mat)
	   (type vec4 vec))
  (with-vec4 vec
	     (macrolet ((row-by-col (row)
				    (declare (type (card8) row))
				    `(+ (the single-float (* (the single-float #{mat ,row 0}) x))
					(the single-float (* (the single-float #{mat ,row 1}) y))
					(the single-float (* (the single-float #{mat ,row 2}) z))
					(the single-float (* (the single-float #{mat ,row 3}) w)))))
	       
	       (setf res (vec4-make (row-by-col 0)
				    (row-by-col 1)
				    (row-by-col 2)
				    (row-by-col 3)))))
  res)
	     

(defun mat44-translate (mat x y z)
  "Multiply the given homogeneous matrix by an implicit translation
matrix. Result stored to same matrix, supposed to be of transform
matrix kind (initialized to identity at the beginning)"
  (declare (type single-float x y z)
	   (type mat44 mat))
  (setf #{mat 0 3} (+ (* x #{mat 0 0})
		      (* y #{mat 0 1})
		      (* z #{mat 0 2})
		      #{mat 0 3})
	#{mat 1 3} (+ (* x #{mat 1 0})
		      (* y #{mat 1 1})
		      (* z #{mat 1 2})
		      #{mat 1 3})
	#{mat 2 3} (+ (* x #{mat 2 0})
		      (* y #{mat 2 1})
		      (* z #{mat 2 2})
		      #{mat 2 3})
	#{mat 3 3} (+ (* x #{mat 3 0})
		      (* y #{mat 3 1})
		      (* z #{mat 3 2})
		      #{mat 3 3}))
  t)
;;possibly: set diag elements to 1 if 0 or unchanged otherwise


(defun mat44-scale (mat x y z)
  "Multiply the given homogeneous matrix by an implicit scaling
matrix. Result stored to given matrix (transform kind)"
  (declare (type single-float x y z)
	   (type mat44 mat))
  (setf #{mat 0 0} (* x #{mat 0 0})
	#{mat 0 1} (* y #{mat 0 1})
	#{mat 0 2} (* z #{mat 0 2})
	#{mat 0 3} (* x #{mat 0 3}) ;; sinon, trans + scal x ne marche pas
	#{mat 1 0} (* x #{mat 1 0})
	#{mat 1 1} (* y #{mat 1 1})
	#{mat 1 2} (* z #{mat 1 2})
	#{mat 1 3} (* y #{mat 1 3}) ;;
	#{mat 2 0} (* x #{mat 2 0})
	#{mat 2 1} (* y #{mat 2 1})
	#{mat 2 2} (* z #{mat 2 2})
	#{mat 2 3} (* z #{mat 2 3}) ;;
	#{mat 3 0} (* x #{mat 3 0})
	#{mat 3 1} (* y #{mat 3 1})
	#{mat 3 2} (* z #{mat 3 2}))
  t)
;;current matrix is a transformation matrix: initially set to identity
;;and progressively modified by transfo matrices.


(defun mat44-rotate (mat angle x y z)
  "Multiply the current homogeneous matrix by an implicit rotation
matrix (rotation of angle around the vector (x,y,z)"
  (declare (type single-float angle x y z)
	   (type mat44 mat))
  (let* ((cosine (cos angle))
	 (sine (sin angle))
	 (cos-diff (- 1 cosine))
	 (xy (* x y))
	 (xz (* x z))
	 (yz (* y z))
	 (xs (* x sine))
	 (ys (* y sine))
	 (zs (* z sine)))
    (declare (type single-float cosine sine cos-diff xy xz yz xs ys zs)
	     (type mat44 mat))
    (mat44-set-null *spare-matrix*)
    (mat44-set *spare-matrix*
	       (+ (* x x cos-diff) cosine)
	       (- (* xy cos-diff) zs)
	       (+ (* xz cos-diff) ys)
	       0.0
	       (+ (* xy cos-diff) zs)
	       (+ (* y y cos-diff) cosine)
	       (- (* yz cos-diff) xs)
	       0.0
	       (- (* xz cos-diff) ys)
	       (+ (* yz cos-diff) xs)
	       (+ (* z z cos-diff) cosine)
	       0.0
	       0.0 0.0 0.0 1.0)
    (mat44-mult mat *spare-matrix* *current-matrix*)))
;;current matrix is a transformation matrix: initially set to identity
;;and progressively modified by transfo matrices.
    


;;pretty-printing facility for displaying mat44 objects
(set-pprint-dispatch
 'mat44
 #'(lambda (s mat)
     (format s "~%~W~%~W~%~W~%~W
GLOS-MATRIX[4*4} OF SINGLE-FLOAT~%"
	     (mat44-row mat 0)
	     (mat44-row mat 1)
	     (mat44-row mat 2)
	     (mat44-row mat 3))))
	     

;;;-----------------------------------------------------------------------------
;;;-----------------------------------------------------------------------------
;; Vectors 
;;;-----------------------------------------------------------------------------
;;;-----------------------------------------------------------------------------


;;;-----------------------------------------------------------------------------
;; Constructors and Accessors
;;;-----------------------------------------------------------------------------

;;We have both 3 and 4 elements vector types
;;Functions or macros names prefixed with "vec-" take normally vec3 arguments. But for commodity,
;;if a vec4 is passed the fourth coordinate w is ignored and the vec4 is treated as a vec3.

(deftype vec3 () '(simple-array single-float (3)))
(deftype vec4 () '(simple-array single-float (4)))


;; Constructors  
(defun vec3-make-null ()
  (make-array 3 :element-type 'single-float))


(defun vec3-make (x y z)
  (declare (type single-float x y z))
  (make-array 3 :element-type 'single-float :initial-contents (list x y z)))


(defun vec4-make-null ()
  (make-array 4 :element-type 'single-float))


(defun vec4-make (x y z w)
  (declare (type single-float x y z w))
  (make-array 4 :element-type 'single-float :initial-contents (list x y z w)))


;;;-----------------------------------------------------------------------------
;; Vector arithmetic operations
;;;-----------------------------------------------------------------------------


(declaim (inline vec3-set vec3-cpy
		 vec4-set vec4-cpy))


(defun vec3-set (v x y z)
  "Set coordinates of a vec3"
  (declare (type vec3 v)
	   (type single-float x y z))
  (setf #{v 0 0} x
	#{v 0 1} y
	#{v 0 2} z)
  t)


(defun vec3-cpy (v1 v2)
  "copy v2 (src) to v1 (dst)"
  (declare (type vec3 v1 v2))
  (setf #{v1 0 0} #{v2 0 0}
	#{v1 0 1} #{v2 0 1}
	#{v1 0 2} #{v2 0 2})
  t)


(defun vec4-set (v x y z w)
  "Set coordinates of a vec4"
  (declare (type vec4 v)
	   (type single-float x y z w))
  (setf #{v 0 0} x
	#{v 0 1} y
	#{v 0 2} z
	#{v 0 3} w)
  t)


(defun vec4-cpy (v1 v2)
  "copy v2 (src) to v1 (dst)"
  (declare (type vec4 v1 v2))
  (setf #{v1 0 0} #{v2 0 0}
	#{v1 0 1} #{v2 0 1}
	#{v1 0 2} #{v2 0 2}
	#{v1 0 3} #{v2 0 3})
  t)


;;;Additional functions for Nurbs

(defun vec4-make-fromvec3 (v)
  "build vec4 instance from vec3 one with homogeneous coordinate to 1.0"
  (vec4-make (vec-x v) (vec-y v) (vec-z v) 1.0))

       
(defun vec3-vector-tovec4 (vec3-array)
  "return a vec4 CL-onedim-array from a vec3 one. w sets to 1.0" 
  (let ((polygon (make-array (length vec3-array) :element-type 'vec4)))
    (map-into polygon #'(lambda (myvec3) (vec4-make-fromvec3 myvec3))
	      vec3-array)))


;;vec4-project defined in VertexBuffer
(defun vec4-project-fromw (v)
  "Project homogeneous space point into cartesian space"
  (let ((w (vec-w v)))
    (if (zerop w)
	+vec3-zero+
	(vec3-make (/ (vec-x v) w)
		   (/ (vec-y v) w)
		   (/ (vec-z v) w)))))


(defun vec4-forcevec3 (v)
  "Brute-force projection of homogeneous point into cartesian space in
just supressing the w coordinate"
  (vec3-make (vec-x v) (vec-y v) (vec-z v)))
  

(defun vec3-round (v)
  "return a vec3 of rounded coordinates"
  (let ((convert (lambda (x) (coerce (round x) 'single-float))))
    (vec3-make (funcall convert (vec-x v))
	       (funcall convert (vec-y v))
	       (funcall convert (vec-z v))))) 


(defun vec4-round (v)
  "return a vec4 of rounded coordinates"
  (let ((convert (lambda (x) (coerce (round x) 'single-float))))
    (vec4-make (funcall convert (vec-x v))
	       (funcall convert (vec-y v))
	       (funcall convert (vec-z v))
	       (funcall convert (vec-w v))))) 
	     

(defun vec3-equal (v1 v2)
  (and (= (vec-x v1) (vec-x v2)) 
       (= (vec-y v1) (vec-y v2))
       (= (vec-z v1) (vec-z v2))))


(defun vec4-equal (v1 v2)
  (and (= (vec-x v1) (vec-x v2)) 
       (= (vec-y v1) (vec-y v2))
       (= (vec-z v1) (vec-z v2))
       (= (vec-w v1) (vec-w v2))))


(defun vec3-add (v1 v2)
  (vec3-make (+ (vec-x v1) (vec-x v2)) 
	     (+ (vec-y v1) (vec-y v2))
	     (+ (vec-z v1) (vec-z v2))))


(defun vec4-add (v1 v2)
  (vec4-make (+ (vec-x v1) (vec-x v2)) 
	     (+ (vec-y v1) (vec-y v2))
	     (+ (vec-z v1) (vec-z v2))
	     (+ (vec-w v1) (vec-w v2))))


(defun vec3-sub (v1 v2)
  (vec3-make (- (vec-x v1) (vec-x v2)) 
	     (- (vec-y v1) (vec-y v2))
	     (- (vec-z v1) (vec-z v2))))


(defun vec4-sub (v1 v2)
  (vec4-make (- (vec-x v1) (vec-x v2)) 
	     (- (vec-y v1) (vec-y v2))
	     (- (vec-z v1) (vec-z v2))
	     (- (vec-w v1) (vec-w v2))))


(defun vec3-abs (v)
  (vec3-make (abs (vec-x v))
	     (abs (vec-y v))
	     (abs (vec-z v))))


(defun vec4-abs (v)
  (vec4-make (abs (vec-x v))
	     (abs (vec-y v))
	     (abs (vec-z v))
	     (abs (vec-w v))))
	     

(defun vec4-addf (v1 v2)
  (setf (vec-x v1) (+ (vec-x v1) (vec-x v2))) 
  (setf (vec-y v1) (+ (vec-y v1) (vec-y v2)))
  (setf (vec-z v1) (+ (vec-z v1) (vec-z v2)))
  (setf (vec-w v1) (+ (vec-w v1) (vec-w v2))))


(defun vec3-mul (v s)
  (declare (type single-float s))
  (vec3-make (* s (vec-x v))
	     (* s (vec-y v))
	     (* s (vec-z v))))


(defun vec3-sqr (v)
  (vec3-make (* (vec-x v) (vec-x v))
	     (* (vec-y v) (vec-y v))
	     (* (vec-z v) (vec-z v))))
	

(defun vec4-sqr (v)
  (vec4-make (* (vec-x v) (vec-x v))
	     (* (vec-y v) (vec-y v))
	     (* (vec-z v) (vec-z v))
	     (* (vec-w v) (vec-w v))))
	

(defun vec4-mul (v s)
  (declare (type single-float s))
  (vec4-make (* s (vec-x v))
	     (* s (vec-y v))
	     (* s (vec-z v))
	     (* s (vec-w v))))


(defun vec3-div (v q)
  (declare (type single-float q))
  (vec3-mul v (/ 1 q)))


(defun vec4-div (v q)
  (declare (type single-float q))
  (vec4-mul v (/ 1 q)))


(defun vec3-maxdeviation (v)
  "return maximum coordinate value in x, y, z"
  (let ((v2 (vec3-abs v)))
    (max (vec-x v2) (vec-y v2) (vec-z v2))))


(defun vec4-maxdeviation (v)
  "return maximum coordinate value in x, y, z, or w"
  (let ((v2 (vec4-abs v)))
    (max (vec-x v2) (vec-y v2) (vec-z v2) (vec-w v2))))
	 

(defun vec-planeangle (v1 v2)
  "compute (radian) angle between given vectors projected to XY plane"
  ;;if angle = pi/2, tan is undefined  
  (- (if (zerop (vec-x v2)) (/ pi 2)
	 (atan (/ (vec-y v2) (vec-x v2))))
     (if (zerop (vec-x v1)) (/ pi 2)
	 (atan (/ (vec-y v1) (vec-x v1))))))
	    

;;these accessors are macros because we want them 'setable'
(defmacro vec-x (v)
  `(aref ,v 0))
(defmacro vec-y (v)
  `(aref ,v 1))
(defmacro vec-z (v)
  `(aref ,v 2))
(defmacro vec-w (v)
  `(aref ,v 3))


(defmacro vec-norm (v)
  "Return the length of the argument vector"
  `(the single-float (sqrt (the single-float (+ (* #{,v 0 0} #{,v 0 0})
						(* #{,v 0 1} #{,v 0 1})
						(* #{,v 0 2} #{,v 0 2}))))))


(defmacro vec-dot-product (v1 v2)
  "Return the dot product of the argument vectors"
  `(the single-float (+ (* #{,v1 0 0} #{,v2 0 0})
			(* #{,v1 0 1} #{,v2 0 1})
			(* #{,v1 0 2} #{,v2 0 2}))))


(defmacro vec-set-cross-product (vtarget v1 v2)
  "vtarget become the result of the cross product between v1 and v2"
  `(vec3-set ,vtarget
	     (- (* #{,v1 0 1} #{,v2 0 2}) (* #{,v1 0 2} #{,v2 0 1}))
	     (- (* #{,v1 0 2} #{,v2 0 0}) (* #{,v1 0 0} #{,v2 0 2}))
	     (- (* #{,v1 0 0} #{,v2 0 1}) (* #{,v1 0 1} #{,v2 0 0}))))


(defmacro vec-set-normalize (v)
  "normalized the vector v and return v"
  (let ((n (gensym)))
    `(let ((,n (vec-norm ,v)))
       (declare (type single-float ,n))
       (if (not (= ,n 0.0))
	   (vec3-set ,v (/ #{,v 0 0} ,n) (/ #{,v 0 1} ,n) (/ #{,v 0 2} ,n))))))


;;IO section


(defun vec3-vector-perimeter (vec3-vector)
  "return the perimeter of given polygon (in ve3 vector form)"
  (let* ((size (length vec3-vector))
	 (prevpoint nil)
	 (distance 0)
	 (mapped (lambda (point)
		   (when prevpoint
		     (incf distance (vec-norm (vec3-sub prevpoint point))))
		   (setf prevpoint point)
		   point)))
    (map-into vec3-vector mapped vec3-vector)
    distance))


(defun vec4-vector-perimeter (vec4-vector)
  "the same accepting vec4 vector input"
  (vec3-vector-perimeter vec4-vector))


(defun vec3-parse-line (string)
  "return a list of vec3 objects from a string text line. If tokens %
3 > 0, remaining tokens are ignored"
  (let ((token-list (str-tokenize string))
	(vec3-objects nil)
	(vec3-coords nil))
    (dolist (token token-list (nreverse vec3-objects))
      (progn (push token vec3-coords)
	     (when (= (length vec3-coords) 3)
	       (let ((vec3-z (str-coerce-sf (car vec3-coords)))
		     (vec3-y (str-coerce-sf (cadr vec3-coords)))
		     (vec3-x (str-coerce-sf (caddr vec3-coords))))
		 (progn (setf vec3-coords nil)
			(push (vec3-make vec3-x vec3-y vec3-z)
			      vec3-objects))))))))


(defun vec3-output-stringize (v)
  "build string for stream output purposes"
  (mkstr (vec-x v) " " (vec-y v) " " (vec-z v)))
  

(defun vec4-output-stringize (v)
  "build string for stream output purposes"
  (mkstr (vec-x v) " " (vec-y v) " " (vec-z v) " " (vec-w v)))
  

;; todo: recognize type to also work with vec4
(defun vec-seq-tofile (vec-seq filename &key (output-dir "./"))
  "write a sequence (vector, list) of vec objects onto given file"
  (let* ((filename (mkstr output-dir filename))
	 (path (make-pathname :name filename)))
    ;;(format t "~%filename: ~A" filename) 
    (with-open-file (stream path :direction :output :if-exists :supersede)
      (map-into vec-seq
		#'(lambda (vec) (progn (write-line (vec3-output-stringize vec) stream)
				       vec))
		vec-seq)))) 


;;pretty-printing facility for displaying vec objects
(set-pprint-dispatch
 'vec3
 #'(lambda (s vec)
     (format s "~%{ ~W; ~W; ~W }
GLOS-VECTOR-3 OF SINGLE-FLOAT~%"
	     (vec-x vec)
	     (vec-y vec)
	     (vec-z vec))))


(set-pprint-dispatch
 'vec4
 #'(lambda (s vec)
     (format s "~%{ ~W; ~W; ~W; ~W }
GLOS-VECTOR-4 OF SINGLE-FLOAT~%"
	     (vec-x vec)
	     (vec-y vec)
	     (vec-z vec)
	     (vec-w vec))))


     
					; x,y,z,w in the body are replaced by #{v 0},#{v 1},...
					; (with-vec4 v (+ x y z w)) -> (+ #{v 0} #{v 1} #{v 2} #{v 3} )
(defmacro with-vec4 (v &rest body)
  (let* ((x `#{,v 0 0})
	 (y `#{,v 0 1})
	 (z `#{,v 0 2})
	 (w `#{,v 0 3})
	 (newbody (substn (list x 'x y 'y z 'z w 'w) body)))
    `(prog ()
	   ,@newbody)))


					; same as with-vec4 with only x,y,z
(defmacro with-vec3 (v &rest body)
  (let* ((x `#{,v 0 0})
	 (y `#{,v 0 1})
	 (z `#{,v 0 2})
	 (newbody (substn (list x 'x y 'y z 'z) body)))
    `(prog ()
	   ,@newbody)))


(defvar +vec3-zero+ (vec3-make-null))
(defvar +vec4-zero+ (vec4-make-null))

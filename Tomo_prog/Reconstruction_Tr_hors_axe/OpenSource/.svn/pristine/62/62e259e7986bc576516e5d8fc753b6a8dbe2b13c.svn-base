;;; ----------------------------------------------------------------------------
;;; $RCSfile: Transforms.lisp,v $
;;; $Author: hubert $
;;; $Date: 2000/05/05 06:13:40 $
;;; $Revision: 1.1.1.1 $
;;; ----------------------------------------------------------------------------

;;GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
;;Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net)


;;;-----------------------------------------------------------------------------
;;; Declaration 
;;;-----------------------------------------------------------------------------


;matrix and vector pools for transform primitives
(defparameter +transforms-matrix-pool+ (make-array 16 :fill-pointer t))
(defparameter +transforms-vector-pool+ (make-array 16 :fill-pointer t))
(defparameter +transforms-matrix-stack-pool+ (make-array 64 :fill-pointer t))


(defparameter *projection-matrix* (mat44-make-null))
(defparameter *modelview-matrix*  (mat44-make-null))
(defparameter *current-matrix* *modelview-matrix*)
(defparameter *current-mode* :MODELVIEW)

;matrix stacks for PushMatrix/PopMatrix
(defparameter *projection-matrix-stack* (make-array 32 :fill-pointer 0))
(defparameter *modelview-matrix-stack* (make-array 32 :fill-pointer 0))


(declaim (type vector
	       +transforms-matrix-pool+
	       +transforms-vector-pool+
	       +transforms-matrix-stack-pool+
	       *projection-matrix-stack*
	       *modelview-matrix-stack*)
	 (type mat44
	       *current-matrix*
	       *projection-matrix*
	       *modelview-matrix*))


;fill pools arrays
(dotimes (i 16)
  (setf (aref +transforms-matrix-pool+ i) (mat44-make-null)
	(aref +transforms-vector-pool+ i) (vec4-make-null)))
(dotimes (i 64)
  (setf (aref +transforms-matrix-stack-pool+ i) (mat44-make-null)))


;;;-----------------------------------------------------------------------------
;;; Matrix stacks stuff
;;;-----------------------------------------------------------------------------


(defun stack-of-current-matrix ()
  "Return the matrix stack of the current mode"
  (if (eq *current-matrix* *modelview-matrix*)
      *modelview-matrix-stack*
    *projection-matrix-stack*))


(defun matrix-of-current-matrix ()
  (if (eq *current-matrix* *modelview-matrix*)
      '*modelview-matrix*
    '*projection-matrix*))


(defmacro MatrixMode (mode)
  "Sets the current matrix mode for future operations on matrices"
  (case mode
    (:MODELVIEW  `(setf *current-matrix* *modelview-matrix*))
    (:PROJECTION `(setf *current-matrix* *projection-matrix*))
    (t (error "Unknown matrix mode: ~S" mode))))
  

(defun PushMatrix ()
  "Push the current-matrix on top of the matrix stack"
  (let ((temp-matrix *current-matrix*)
	(matrix-of-current-mode (matrix-of-current-matrix)))
    (declare (type mat44 temp-matrix)
	     (type symbol matrix-of-current-mode))
    (vector-push *current-matrix* (stack-of-current-matrix))
    (setf *current-matrix* (vector-pop +transforms-matrix-stack-pool+))
    (mat44-copy *current-matrix* temp-matrix)
    (set matrix-of-current-mode *current-matrix*)
    t))


(defun PopMatrix ()
  "Set the current matrix to the matrix on top of the matrix stack"
  (let ((matrix-of-current-mode (matrix-of-current-matrix)))
    (declare (type symbol matrix-of-current-mode))
    (vector-push *current-matrix* +transforms-matrix-stack-pool+)
    (setf *current-matrix* (vector-pop (stack-of-current-matrix)))
    (set matrix-of-current-mode *current-matrix*)
    t))

;;;-----------------------------------------------------------------------------
;;; Matrix Transformations
;;;-----------------------------------------------------------------------------


; on ne peut pas utiliser le reader macro ici a cause de l'ordre de la macroexpansion:
; la reader macro est toujours expansee *avant* le reste de la macro. donc on utilise mat44-elem.
(defmacro row-by-col (mat-dest mat1 mat2 i j)
  (declare (type fixnum i j))
  `(setf (mat44-elem ,mat-dest ,i ,j)
	 (the single-float (+ (the single-float (* (mat44-elem ,mat1 ,i 0) (mat44-elem ,mat2 0 ,j)))
			      (the single-float (* (mat44-elem ,mat1 ,i 1) (mat44-elem ,mat2 1 ,j)))
			      (the single-float (* (mat44-elem ,mat1 ,i 2) (mat44-elem ,mat2 2 ,j)))
			      (the single-float (* (mat44-elem ,mat1 ,i 3) (mat44-elem ,mat2 3 ,j)))))))
			      

(defmacro calc-line (mat-dest mat1 mat2 i)
  (declare (type fixnum i))
  `(progn (row-by-col ,mat-dest ,mat1 ,mat2 ,i 0)
	  (row-by-col ,mat-dest ,mat1 ,mat2 ,i 1)
	  (row-by-col ,mat-dest ,mat1 ,mat2 ,i 2)
	  (row-by-col ,mat-dest ,mat1 ,mat2 ,i 3)))


(defun MultMatrix (mat)
  "Multiply the current matrix by the argument matrix and store the result into the current matrix"
   (let ((m (vector-pop +transforms-matrix-pool+)))
     (declare (type mat44 mat m))
     (mat44-copy m *current-matrix*)
     (calc-line *current-matrix* m mat 0)
     (calc-line *current-matrix* m mat 1)
     (calc-line *current-matrix* m mat 2)
     (calc-line *current-matrix* m mat 3)
     (vector-push m +transforms-matrix-pool+)))
				  
			  
(defmacro LoadMatrix (m)
  "Set the current matrix to an arbitrary matrix"
  `(mat44-copy *current-matrix* ,m))


(defun LoadIdentity ()
  "Set the current matrix to the identity matrix"
  (mat44-set *current-matrix*
	     1.0 0.0 0.0 0.0
	     0.0 1.0 0.0 0.0
	     0.0 0.0 1.0 0.0
	     0.0 0.0 0.0 1.0)
  t)


(defun Translate (x y z)
  (declare (type single-float x y z))
  "Multiply the current matrix by a translation matrix"
  (setf #{*current-matrix* 0 3} (+ (* x #{*current-matrix* 0 0})
				   (* y #{*current-matrix* 0 1})
				   (* z #{*current-matrix* 0 2})
				   #{*current-matrix* 0 3})
	#{*current-matrix* 1 3} (+ (* x #{*current-matrix* 1 0})
				   (* y #{*current-matrix* 1 1})
				   (* z #{*current-matrix* 1 2})
				   #{*current-matrix* 1 3})
	#{*current-matrix* 2 3} (+ (* x #{*current-matrix* 2 0})
				   (* y #{*current-matrix* 2 1})
				   (* z #{*current-matrix* 2 2})
				   #{*current-matrix* 2 3})
	#{*current-matrix* 3 3} (+ (* x #{*current-matrix* 3 0})
				   (* y #{*current-matrix* 3 1})
				   (* z #{*current-matrix* 3 2})
				   #{*current-matrix* 3 3}))
  t)


(defun Scale (x y z)
  "Multiply the current matrix by a scaling matrix"
  (declare (type single-float x y z))
  (setf #{*current-matrix* 0 0} (* x #{*current-matrix* 0 0})
	#{*current-matrix* 0 1} (* y #{*current-matrix* 0 1})
	#{*current-matrix* 0 2} (* z #{*current-matrix* 0 2})
	#{*current-matrix* 1 0} (* x #{*current-matrix* 1 0})
	#{*current-matrix* 1 1} (* y #{*current-matrix* 1 1})
	#{*current-matrix* 1 2} (* z #{*current-matrix* 1 2})
	#{*current-matrix* 2 0} (* x #{*current-matrix* 2 0})
	#{*current-matrix* 2 1} (* y #{*current-matrix* 2 1})
	#{*current-matrix* 2 2} (* z #{*current-matrix* 2 2})
	#{*current-matrix* 3 0} (* x #{*current-matrix* 3 0})
	#{*current-matrix* 3 1} (* y #{*current-matrix* 3 1})
	#{*current-matrix* 3 2} (* z #{*current-matrix* 3 2}))
  t)


(defun Rotate (angle x y z)
  "Multiply the current matrix by a rotation matrix (rotation of angle around the vector (x,y,z)"
  (declare (type single-float angle x y z))
  (let* ((m (vector-pop +transforms-matrix-pool+))
	  (cosine (cos angle))
	  (sine (sin angle))
	  (cos-diff (- 1 cosine))
	  (xy (* x y))
	  (xz (* x z))
	  (yz (* y z))
	  (xs (* x sine))
	  (ys (* y sine))
	  (zs (* z sine)))
    (declare (type single-float cosine sine cos-diff xy xz yz xs ys zs)
	     (type mat44 m))
    (mat44-set m
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
    (MultMatrix m)
    (vector-push m +transforms-matrix-pool+)))


(defun Frustum (left right bottom top near far)
  "Set the frustum"
  (declare (type single-float left right bottom top near far))          
  (let*  ((mrl (- right left))
	 (mtb (- top bottom))
	 (mfn (- far near))
	 (m00 (/ (* 2.0 near) mrl))
	 (m02 (/ (+ right left) mrl))
	 (m11 (/ (* 2.0 near) mtb))
	 (m12 (/ (+ top bottom) mtb))
	 (m22 (/ (- (+ far near)) mfn))
	 (m23 (/ (- (* 2.0 far near)) mfn))
	 (frustum-matrix (vector-pop +transforms-matrix-pool+)))
    (declare (type single-float mrl mtb mfn m00 m02 m11 m12 m22 m23)
	     (type mat44 frustum-matrix))
    (setf *near* (the card32 (truncate (the single-float near))))
    (setf *far* (the card32 (truncate (the single-float far))))
    (mat44-set frustum-matrix
	       m00 0.0 m02 0.0
	       0.0 m11 m12 0.0
	       0.0 0.0 m22 m23
	       0.0 0.0 -1.0 0.0)
    (MultMatrix frustum-matrix)
    (vector-push frustum-matrix +transforms-matrix-pool+)))


(defun Perspective (fovy aspect znear zfar)
  "Another way to set the frustum using the fov"
  (declare (type single-float fovy aspect znear zfar))
  (let* ((ymax (* znear (tan (* fovy (/ +PI+ 360.0)))))
	 (ymin (- ymax)))
    (declare (type single-float ymin ymax))
    (Frustum (* ymin aspect) (* ymax aspect) ymin ymax znear zfar)))


(defun LookAt (eyeX eyeY eyeZ centerX centerY centerZ upX upY upZ)
  "Set the viewer to (eyeX,eyeY,eyeZ) looking at (centerX,centerY,centerZ) with up vector (upX,upY,upZ)"
  (declare (type single-float eyeX eyeY eyeZ centerX centerY centerZ upX upY upZ))
  (let ((X (vector-pop +transforms-vector-pool+))
	(Y (vector-pop +transforms-vector-pool+))
	(Z (vector-pop +transforms-vector-pool+))
	(m (vector-pop +transforms-matrix-pool+)))
    (declare (type vec4 X Y Z)
	     (type mat44 m))
    (vec4-set Z (- eyeX centerX) (- eyeY centerY) (- eyeZ CenterZ) 1.0)
    (vec-set-normalize Z)
    (vec4-set Y upX upY upZ 1.0)
    (vec-set-cross-product X Y Z)
    (vec-set-cross-product Y Z X)
    (vec-set-normalize X)
    (vec-set-normalize Y)
    (mat44-set m 
	       #{X 0} #{Y 0} #{Z 0} 0.0
	       #{X 1} #{Y 1} #{Z 1} 0.0
	       #{X 2} #{Y 2} #{Z 2} 0.0
	       0.0    0.0    0.0    1.0)
    (MultMatrix m)
    (Translate (- eyeX) (- eyeY) (- eyeZ))
    (vector-push X +transforms-vector-pool+)
    (vector-push Y +transforms-vector-pool+)
    (vector-push Z +transforms-vector-pool+)
    (vector-push m +transforms-matrix-pool+)))

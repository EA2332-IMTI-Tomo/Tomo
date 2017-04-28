;;; ----------------------------------------------------------------------------
;;; $RCSfile: DemoASC.lisp,v $
;;; $Author: hubert $
;;; $Date: 2000/05/05 06:13:41 $
;;; $Revision: 1.1.1.1 $
;;; ----------------------------------------------------------------------------

;;GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
;;Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net)


;;must be done without compilation
;;(load "./autorun-clx.lisp")


(defmacro with-display-window (body &key (host "localhost")
				    (minimized-title "GLOS display")
				    (window-title "GLOS display window")
				    (font "fixed")
				    (scr-width 640)
				    (scr-height 400))
  "define an environment opening an X Window of &key attribute,
closable by single lmouse click. Given body code will be processed
inside this environment, where the function (put-point (x y)) is
valid (x and y integers not necessarily bounded)."
  `(let ((display nil)
	 (abort t))
     (unwind-protect
	 (progn 
	   (setf display (open-display ,host))
	   (multiple-value-prog1
	    (let* ((screen (display-default-screen display))
		   (black (screen-black-pixel screen))
		   (white (screen-white-pixel screen))
		   (font (open-font display ,font))
		   (width ,scr-width)
		   (height ,scr-height)
		   (x (truncate (- (screen-width screen) width) 2))
		   (y (truncate (- (screen-height screen) height) 2))
		   (window (create-window :parent (screen-root screen)
					  :x x :y y :width width :height height
					  :background black
					  :border white
					  :border-width 1
					  :colormap (screen-default-colormap screen)
					  :bit-gravity :center
					  :event-mask '(:exposure :button-press)))
		   (gcontext (create-gcontext :drawable window
					      :background black
					      :foreground white
					      :font font)))
	      ;; Set window manager hints
	      (set-wm-properties window
				 :name ,window-title
				 :icon-name ,minimized-title
				 :resource-name ,minimized-title
				 :resource-class 'Test-GLOS
				 :command (list 'Test-GLOS ,host)
				 :x x :y y :width width :height height
				 :min-width width :min-height height
				 :input :off :initial-state :normal)
	     
	      (map-window window) 
	      ;; Handle events
	      (event-case (display :discard-p t :force-output-p t)
			  (exposure;; Come here on exposure events
			   (window count)

			   (flet ((put-point (x y)
					     (draw-point window gcontext x y)))
			     
			     (progn
			       ;; insert display code here!!!
			       ,body
			       nil	; prevents from returning immediately 
			       ;; end of display code
			       )))
	      
			  (button-press () t))
	      (setf abort nil))))
       ;; Ensure display is closed when done
       (when display
	 (close-display display :abort abort)))))






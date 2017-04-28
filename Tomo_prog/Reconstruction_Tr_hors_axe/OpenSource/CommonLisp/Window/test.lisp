(with-display-window
 (labels ((vec3-window-plot (v)
			    (let ((x (round (vec-x v)))
				  (y (round (vec-y v))))
			      (put-point x y)))
	  (vec4-window-plot (v) (vec3-window-plot v)))
   (progn
     (put-point 12 23)
     (put-point 122 231)))
 :host "hydra")

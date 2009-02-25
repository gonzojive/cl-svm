(cl:in-package #:cl-user)

(defpackage #:cl-svm-system
  (:use #:cl #:asdf))

(in-package :cl-svm-system)

(defsystem :cl-svm
  :name "cl-svm"
  :description ""
  :version "02.09"
  :author "Red Daly <reddaly@gmail.com>"
  :license "Communist license"
  :components ((:static-file "cl-svm.asd")
	       (:module
		:src
		:components ((:file "package")
			     (:file "svm" :depends-on ("package")))))
  :depends-on ())
			
			     
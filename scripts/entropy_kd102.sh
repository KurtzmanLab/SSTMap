#!/bin/bash

         for i in $( ls cluster.*pdb );
         do
         ./kdhsa102 -e $i -i ../organizewaters_4/$i
	 #echo $i
         done



 

/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */


package org.biojava.bio.seq;

/**
 * Title:        Frame.<p>
 * Description:  An interface that can be implemented by classes that contain
 *               frame information.<p>
 * Copyright:    Copyright (c) 2001.<p>
 * @author Mark Schreiber
 * @since 1.2
 * @version 1.0
 */

public interface Frame {
  public static int FRAME_0 = 0;
  public static int FRAME_1 = 1;
  public static int FRAME_2 = 2;

  /**
   * A method to get the frame information of the implementing object
   * @return an integer from 0 to 2 representing the frame.
   */
  public int getFrame();
}

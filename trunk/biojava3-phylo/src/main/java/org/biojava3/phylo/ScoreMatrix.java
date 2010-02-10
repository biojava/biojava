/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.phylo;
/*
 * This source file is derived from jalview.schemes.ScoreMatrix to provide
 * tree construction in the biojava package to minimize the number of external
 * JalView classes required.
 *
 * Jalview - A Sequence Alignment Editor and Viewer (Version 2.4)
 * Copyright (C) 2008 AM Waterhouse, J Procter, G Barton, M Clamp, S Searle
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */


public class ScoreMatrix
{
  String name;

  /**
   * reference to integer score matrix
   */
  int[][] matrix;

  /**
   * 0 for Protein Score matrix. 1 for dna score matrix
   */
  int type;

  ScoreMatrix(String name, int[][] matrix, int type)
  {
    this.matrix = matrix;
    this.type = type;
  }

  public boolean isDNA()
  {
    return type == 1;
  }

  public boolean isProtein()
  {
    return type == 0;
  }

  public int[][] getMatrix()
  {
    return matrix;
  }

  /**
   *
   * @param A1
   * @param A2
   * @return score for substituting first char in A1 with first char in A2
   */
  public int getPairwiseScore(String A1, String A2)
  {
    return getPairwiseScore(A1.charAt(0), A2.charAt(0));
  }

  


  public int getPairwiseScore(char c, char d)
  {
    int pog = 0;

    try
    {
      int a = (type == 0) ? ResidueProperties.aaIndex[c]
              : ResidueProperties.nucleotideIndex[c];
      int b = (type == 0) ? ResidueProperties.aaIndex[d]
              : ResidueProperties.nucleotideIndex[d];

      pog = matrix[a][b];
    } catch (Exception e)
    {
      // System.out.println("Unknown residue in " + A1 + " " + A2);
    }

    return pog;
  }

}


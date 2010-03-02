/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.exceptions;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class HeaderParseException extends Error{

  private static final long serialVersionUID = -8356845980320906455L;

    public HeaderParseException(String error){
        super(error);
    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.exceptions;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FileAccessError extends Error{

    public FileAccessError(String error){
        super(error);
    }
}

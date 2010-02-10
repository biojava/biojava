/*
 *                  BioJava development code
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
 * Created on Jan 3, 2008
 * 
 */


package org.biojava.dasobert.das;

import java.util.*;

/**
 * Status codes for DAS responses.
 *
 * @author Thomas Down
 * @since 1.1
 */

public class DASStatus {
    public static final int STATUS_OKAY = 200;
    public static final int STATUS_OK = 200;

    public static final int STATUS_BAD_COMMAND = 400;
    public static final int STATUS_BAD_DATASOURCE = 401;
    public static final int STATUS_BAD_COMMAND_ARGUMENTS = 402;
    public static final int STATUS_BAD_REFERENCE = 403;
    public static final int STATUS_BAD_STYLESHEET = 404;
    public static final int STATUS_BAD_COORDS = 405;

    public static final int STATUS_SERVER_ERROR = 500;
    public static final int STATUS_UNSUPPORTED_FEATURE = 501;
    
    private static final Map<Integer,String> errorMessages;

    static {
        errorMessages = new HashMap<Integer,String>();
        errorMessages.put(new Integer(STATUS_BAD_COMMAND), "Bad command");
        errorMessages.put(new Integer(STATUS_BAD_DATASOURCE), "Bad datasource");
        errorMessages.put(new Integer(STATUS_BAD_COMMAND_ARGUMENTS), "Bad command arguments");
        errorMessages.put(new Integer(STATUS_BAD_REFERENCE), "Bad reference");
        errorMessages.put(new Integer(STATUS_BAD_STYLESHEET), "Bad stylesheet");
        errorMessages.put(new Integer(STATUS_BAD_COORDS), "Bad coordinates");
        errorMessages.put(new Integer(STATUS_SERVER_ERROR), "Server error");
        errorMessages.put(new Integer(STATUS_UNSUPPORTED_FEATURE), "Unimplemented feature");
    }
    
    public static String UNKNOWN = "Unknown error";
    
    public static String getErrorDescription(int code) {
        String desc = (String) errorMessages.get(new Integer(code));
        if (desc == null) {
            desc = UNKNOWN;
        }
        return desc;
    }
}

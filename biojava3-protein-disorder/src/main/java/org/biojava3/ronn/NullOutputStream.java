/* 
 * @(#)NullOutputStream.java	1.0 June 2010
 * 
 * Copyright (c) 2010 Peter Troshin
 * 
 * JRONN version: 3.1     
 *  
 *  This library is free software; you can redistribute it and/or modify it under the terms of the
 *  Apache License version 2 as published by the Apache Software Foundation
 * 
 *  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 *  even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache 
 *  License for more details.
 * 
 *  A copy of the license is in apache_license.txt. It is also available here:
 * see: http://www.apache.org/licenses/LICENSE-2.0.txt
 * 
 * Any republication or derived work distributed in source code form
 * must include this copyright and license notice.
 */
package org.biojava3.ronn;

import java.io.IOException;
import java.io.OutputStream;

/**
 * The stream that void its input
 * 
 * @author Petr Troshin
 * @version 1.0
 * @since 3.0.2
 */
public final class NullOutputStream extends OutputStream {

    @Override
    public void write(final int b) throws IOException {
	// this methods does nothing.
	// This is an intention
    }

}

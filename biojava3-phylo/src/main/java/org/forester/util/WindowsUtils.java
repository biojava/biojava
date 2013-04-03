// $Id: WindowsUtils.java,v 1.2 2009/10/26 23:29:40 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// From: http://www.rgagnon.com/javadetails/java-0652.html
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.util;

import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;

public class WindowsUtils {

    private static final String REGQUERY_UTIL      = "reg query ";
    private static final String REGSTR_TOKEN       = "REG_SZ";
    private static final String DESKTOP_FOLDER_CMD = REGQUERY_UTIL
                                                           + "\"HKCU\\Software\\Microsoft\\Windows\\CurrentVersion\\"
                                                           + "Explorer\\Shell Folders\" /v DESKTOP";

    private WindowsUtils() {
    }

    public static String getCurrentUserDesktopPath() {
        try {
            final Process process = Runtime.getRuntime().exec( DESKTOP_FOLDER_CMD );
            final StreamReader reader = new StreamReader( process.getInputStream() );
            reader.start();
            process.waitFor();
            reader.join();
            final String result = reader.getResult();
            final int p = result.indexOf( REGSTR_TOKEN );
            if ( p == -1 ) {
                return null;
            }
            return result.substring( p + REGSTR_TOKEN.length() ).trim();
        }
        catch ( final Exception e ) {
            return null;
        }
    }

    static class StreamReader extends Thread {

        private final InputStream  is;
        private final StringWriter sw;

        StreamReader( final InputStream is ) {
            this.is = is;
            sw = new StringWriter();
        }

        String getResult() {
            return sw.toString();
        }

        @Override
        public void run() {
            try {
                int c;
                while ( ( c = is.read() ) != -1 ) {
                    sw.write( c );
                }
            }
            catch ( final IOException e ) {
                // Do nothing
            }
        }
    }
}

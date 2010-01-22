// $Id: Printer.java,v 1.1 2009/02/22 00:56:05 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2009-2010 Christian M. Zmasek
// Copyright (C) 2009-2010 Burnham Institute for Medical Research
// All rights reserved
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

package org.forester.archaeopteryx;

import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;

import org.forester.util.ForesterUtil;

final class Printer {

    private Printer() {
        // Hidden constructor.
    }

    /**
     * Returns null if printing has been aborted by the user,
     * a String otherwise -- if a printer name was obtained this String is 
     * the printer name, an empty String otherwise.
     * 
     * @param tree_panel
     * @param job_name
     * @return
     * @throws PrinterException
     */
    static String print( final TreePanel tree_panel, final String job_name ) throws PrinterException {
        if ( ( tree_panel == null ) || ( tree_panel.getPhylogeny() == null ) ) {
            throw new IllegalArgumentException( "attempt to print null" );
        }
        if ( ForesterUtil.isEmpty( job_name ) ) {
            throw new IllegalArgumentException( "attempt use null or empty print job name" );
        }
        final PrinterJob printer_job = PrinterJob.getPrinterJob();
        if ( printer_job != null ) {
            printer_job.setJobName( job_name );
            printer_job.setPrintable( tree_panel );
            final boolean ok = printer_job.printDialog();
            if ( ok ) {
                printer_job.print();
                final String print_service_name = printer_job.getPrintService().getName();
                if ( !ForesterUtil.isEmpty( print_service_name ) ) {
                    return print_service_name;
                }
                return "";
            }
            else {
                return null;
            }
        }
        else {
            throw new PrinterException( "failed to access printer job" );
        }
    }
}
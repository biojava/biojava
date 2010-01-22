// $Id: PdfWriter.java,v 1.3 2009/02/23 18:59:17 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
// Copyright (C) 2003-2007 Ethalinda K.S. Cannon
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
//removed to deal with Biojava Maven integration
//import gnu.jpdf.PDFGraphics;
//import gnu.jpdf.PDFJob;

import java.awt.Graphics;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import org.forester.phylogeny.Phylogeny;

/*
 * 
 * This uses gnujpdf. Version 1.7.0 seems to produce pdf-files with errors when
 * using Linux ggv or acroread (7.2).
 * 
 * See: http://sourceforge.net/projects/gnujpdf/ and
 * http://gnujpdf.sourceforge.net/
 */
final class PdfWriter {

    private PdfWriter() {
        // Empty constructor.
    }

    public static String writePhylogenyToPdf( final String file_name, final TreePanel tree_panel, final Options options )
            throws IOException {
        return writePhylogenyToPdf( file_name, tree_panel, options, 0, 0 );
    }

    public static String writePhylogenyToPdf( final String file_name,
                                              final TreePanel tree_panel,
                                              final Options options,
                                              final int width,
                                              final int height ) throws IOException {
        final Phylogeny phylogeny = tree_panel.getPhylogeny();
        if ( ( phylogeny == null ) || phylogeny.isEmpty() ) {
            return "";
        }
        final File file = new File( file_name );
        if ( file.isDirectory() ) {
            throw new IllegalArgumentException( "[" + file_name + "] is a directory" );
        }
        final FileOutputStream os = new FileOutputStream( file );
        
//removed to deal with BioJava integration
//        final PDFJob job = new PDFJob( os );
//        if ( job == null ) {
//            throw new IOException( "failed to obtain PDF job" );
//        }
//        final Graphics pdf_graphics = job.getGraphics();
//        if ( pdf_graphics == null ) {
//            throw new IOException( "failed to obtain PDF graphics" );
//        }
//        ( ( PDFGraphics ) pdf_graphics ).setLineWidth( options.getPrintLineWidth() );
//        if ( ( width > 0 ) && ( height > 0 ) ) {
//            try {
//                job.getPageDimension().setSize( width, height );
//            }
//            catch ( final Exception e ) {
//                throw new IOException( "failed to change pdf output size: " + e );
//            }
//        }
//        tree_panel.paintPhylogeny( pdf_graphics, true, false, 0, height, 0, 0 );
//        pdf_graphics.dispose();
//        job.end();
        os.flush();
        os.close();
        String msg = file.toString();
        if ( ( width > 0 ) && ( height > 0 ) ) {
            msg += " [size: " + width + ", " + height + "]";
        }
        return msg;
    }
}

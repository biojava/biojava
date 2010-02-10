// $Id: PdfExporter.java,v 1.6 2009/11/15 08:41:19 cmzmasek Exp $
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

import java.awt.Graphics2D;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import org.forester.phylogeny.Phylogeny;
//removed to deal with biojava integration
//import com.lowagie.text.Document;
//import com.lowagie.text.DocumentException;
//import com.lowagie.text.FontFactory;
//import com.lowagie.text.Rectangle;
//import com.lowagie.text.pdf.DefaultFontMapper;
//import com.lowagie.text.pdf.PdfContentByte;
//import com.lowagie.text.pdf.PdfWriter;

/*
 * 
 * This uses iText.
 * 
 * See: http://www.lowagie.com/iText/
 * 
 * Current version: iText-2.1.7
 */
final class PdfExporter {

    private static final int HEIGHT_LIMIT = 100;
    private static final int WIDTH_LIMIT  = 60;

    private PdfExporter() {
        // Empty constructor.
    }

    static String writePhylogenyToPdf( final String file_name, final TreePanel tree_panel, int width, int height )
            throws IOException {
        if ( height < HEIGHT_LIMIT ) {
            height = HEIGHT_LIMIT;
        }
        if ( width < WIDTH_LIMIT ) {
            width = WIDTH_LIMIT;
        }
        final Phylogeny phylogeny = tree_panel.getPhylogeny();
        if ( ( phylogeny == null ) || phylogeny.isEmpty() ) {
            return "";
        }
        if ( tree_panel.getMainPanel().getTreeFontSet().getSmallFont().getSize() < 1 ) {
            throw new IOException( "fonts are too small for PDF export" );
        }
        final File file = new File( file_name );
        if ( file.isDirectory() ) {
            throw new IOException( "[" + file_name + "] is a directory" );
        }
/*        final Document document = new Document();
//        document.setPageSize( new Rectangle( width, height ) );
//        document.setMargins( WIDTH_LIMIT / 2, WIDTH_LIMIT / 2, HEIGHT_LIMIT / 2, HEIGHT_LIMIT / 2 );
        PdfWriter writer = null;
        try {
            writer = PdfWriter.getInstance( document, new FileOutputStream( file_name ) );
        }
        catch ( final DocumentException e ) {
            throw new IOException( e );
        }
        document.open();
        final DefaultFontMapper mapper = new DefaultFontMapper();
        FontFactory.registerDirectories();
        if ( Util.isWindows() ) {
            mapper.insertDirectory( "C:\\WINDOWS\\Fonts\\" );
        }
        else if ( Util.isMac() ) {
            mapper.insertDirectory( "/Library/Fonts/" );
            mapper.insertDirectory( "/System/Library/Fonts/" );
        }
        else {
            mapper.insertDirectory( "/usr/X/lib/X11/fonts/TrueType/" );
            mapper.insertDirectory( "/usr/X/lib/X11/fonts/Type1/" );
            mapper.insertDirectory( "/usr/share/fonts/default/TrueType/" );
            mapper.insertDirectory( "/usr/share/fonts/default/Type1/" );
        }
        final PdfContentByte cb = writer.getDirectContent();
        final Graphics2D g2 = cb.createGraphics( width, height, mapper );
        try {
            tree_panel.paintPhylogeny( g2, true, false, width - WIDTH_LIMIT, height - HEIGHT_LIMIT, 0, 0 );
        }
        catch ( final Exception e ) {
            Util.unexpectedException( e );
        }
        finally {
            try {
                g2.dispose();
                document.close();
            }
            catch ( final Exception e ) {
                //Do nothing.
            }
        }
 * 
 */
        String msg = file.toString();
        if ( ( width > 0 ) && ( height > 0 ) ) {
            msg += " [size: " + width + ", " + height + "]";
        }
        return msg;
    }
}

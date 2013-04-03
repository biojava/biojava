// $Id: BranchColor.java,v 1.14 2009/01/13 19:49:29 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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

package org.forester.phylogeny.data;

import java.awt.Color;
import java.io.IOException;
import java.io.Writer;

import org.forester.io.parsers.nhx.NHXtags;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.util.ForesterUtil;

public class BranchColor implements PhylogenyData {

    private Color _color;

    public BranchColor() {
        _color = null;
    }

    public BranchColor( final Color color ) {
        _color = color;
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( getValue().toString() );
    }

    @Override
    public StringBuffer asText() {
        return new StringBuffer( getValue().toString() );
    }

    @Override
    /**
     * Not a deep copy.
     * 
     */
    public PhylogenyData copy() {
        final BranchColor bc = new BranchColor();
        bc.setValue( getValue() );
        return bc;
    }

    public Color getValue() {
        return _color;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        return getValue().equals( ( ( BranchColor ) data ).getValue() );
    }

    public void setValue( final Color color ) {
        _color = color;
    }

    @Override
    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        sb.append( NHXtags.COLOR );
        sb.append( getValue().getRed() );
        sb.append( "." );
        sb.append( getValue().getGreen() );
        sb.append( "." );
        sb.append( getValue().getBlue() );
        return sb;
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.COLOR );
        PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.COLOR_RED, getValue().getRed() + "", indentation );
        PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.COLOR_GREEN, getValue().getGreen() + "", indentation );
        PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.COLOR_BLUE, getValue().getBlue() + "", indentation );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.COLOR );
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}

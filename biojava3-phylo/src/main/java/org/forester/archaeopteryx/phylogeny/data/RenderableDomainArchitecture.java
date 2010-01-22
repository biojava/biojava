// $Id: RenderableDomainArchitecture.java,v 1.2 2009/01/13 19:24:57 cmzmasek Exp
// $
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

package org.forester.archaeopteryx.phylogeny.data;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.io.IOException;
import java.io.Writer;
import java.util.Map;
import java.util.SortedMap;

import org.forester.archaeopteryx.Constants;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.util.ForesterUtil;

public final class RenderableDomainArchitecture extends DomainArchitecture implements RenderablePhylogenyData {

    static private Map<String, Color> Domain_colors;
    final static private int          BRIGHTEN_COLOR_BY             = 200;
    final static private int          E_VALUE_THRESHOLD_EXP_DEFAULT = 0;
    private static int                _Next_default_domain_color    = 0;
    private final static String[]     DEFAULT_DOMAINS_COLORS        = { "0xFF0000", "0x0000FF", "0xAAAA00", "0xFF00FF",
            "0x00FFFF", "0x800000", "0x000080", "0x808000", "0x800080", "0x008080", "0xE1B694" };
    private int                       _e_value_threshold_exp        = RenderableDomainArchitecture.E_VALUE_THRESHOLD_EXP_DEFAULT;
    private boolean                   _display_domain_names         = true;
    private double                    _rendering_factor_width       = 1.0;
    private double                    _rendering_height             = 0;
    private final DomainArchitecture  _domain_structure;
    private final Rectangle2D         _rectangle                    = new Rectangle2D.Float();

    public RenderableDomainArchitecture( final DomainArchitecture domain_structure ) {
        _domain_structure = domain_structure;
    }

    @Override
    public StringBuffer asSimpleText() {
        return _domain_structure.asSimpleText();
    }

    @Override
    public StringBuffer asText() {
        return _domain_structure.asText();
    }

    @Override
    public PhylogenyData copy() {
        return _domain_structure.copy();
    }

    private final void drawDomain( final double x,
                                   final double y,
                                   final double width,
                                   final double heigth,
                                   final String name,
                                   final Graphics2D g,
                                   final boolean to_pdf ) {
        final double h2 = heigth / 2.0;
        final Color color_one = getColorOne( name );
        final Color color_two = getColorTwo( color_one );
        double step = 1;
        if ( to_pdf ) {
            step = 0.1;
        }
        for( double i = 0; i < heigth; i += step ) {
            g.setColor( org.forester.util.ForesterUtil
                    .calcColor( i >= h2 ? heigth - i : i, 0, h2, color_one, color_two ) );
            _rectangle.setFrame( x, i + y, width, step );
            g.fill( _rectangle );
        }
    }

    private Color getColorOne( final String name ) {
        Color c = Constants.DOMAIN_STRUCTURE_COLOR_1;
        if ( RenderableDomainArchitecture.Domain_colors != null ) {
            c = RenderableDomainArchitecture.Domain_colors.get( name );
            if ( c == null ) {
                if ( RenderableDomainArchitecture._Next_default_domain_color < RenderableDomainArchitecture.DEFAULT_DOMAINS_COLORS.length ) {
                    c = Color
                            .decode( RenderableDomainArchitecture.DEFAULT_DOMAINS_COLORS[ RenderableDomainArchitecture._Next_default_domain_color++ ] );
                    RenderableDomainArchitecture.Domain_colors.put( name, c );
                }
                else {
                    c = Constants.DOMAIN_STRUCTURE_COLOR_1;
                }
            }
        }
        return c;
    }

    private Color getColorTwo( final Color color_one ) {
        final int red = color_one.getRed() + RenderableDomainArchitecture.BRIGHTEN_COLOR_BY;
        final int green = color_one.getGreen() + RenderableDomainArchitecture.BRIGHTEN_COLOR_BY;
        final int blue = color_one.getBlue() + RenderableDomainArchitecture.BRIGHTEN_COLOR_BY;
        return new Color( red > 255 ? 255 : red, green > 255 ? 255 : green, blue > 255 ? 255 : blue );
    }

    @Override
    public ProteinDomain getDomain( final int i ) {
        return _domain_structure.getDomain( i );
    }

    @Override
    public SortedMap<Double, ProteinDomain> getDomains() {
        return _domain_structure.getDomains();
    }

    @Override
    public int getNumberOfDomains() {
        return _domain_structure.getNumberOfDomains();
    }

    public Dimension getOriginalSize() {
        return new Dimension( _domain_structure.getTotalLength(), ForesterUtil.roundToInt( _rendering_height ) );
    }

    public Object getParameter() {
        return new Integer( _e_value_threshold_exp );
    }

    public double getRenderingFactorWidth() {
        return _rendering_factor_width;
    }

    public Dimension getRenderingSize() {
        return new Dimension( ForesterUtil.roundToInt( _domain_structure.getTotalLength() * _rendering_factor_width ),
                              ForesterUtil.roundToInt( _rendering_height ) );
    }

    @Override
    public int getTotalLength() {
        return _domain_structure.getTotalLength();
    }

    public boolean isDisplayDomainNames() {
        return _display_domain_names;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        return _domain_structure.isEqual( data );
    }

    public void render( final double x1,
                        final double y1,
                        final Graphics2D g,
                        final TreePanel tree_panel,
                        final boolean to_pdf ) {
        final double f = getRenderingFactorWidth();
        final double y = y1 + ( _rendering_height / 2 );
        final double start = x1 + 20.0;
        g.setColor( Constants.DOMAIN_STRUCTURE_FONT_COLOR );
        _rectangle.setFrame( start, y - 0.5, _domain_structure.getTotalLength() * f, 1 );
        g.fill( _rectangle );
        for( int i = 0; i < _domain_structure.getDomains().size(); ++i ) {
            final ProteinDomain d = _domain_structure.getDomain( i );
            if ( d.getConfidence() <= Math.pow( 10, _e_value_threshold_exp ) ) {
                final double xa = start + d.getFrom() * f;
                final double xb = xa + d.getLength() * f;
                if ( isDisplayDomainNames() ) {
                    g.setFont( tree_panel.getMainPanel().getTreeFontSet().getSmallFont() );
                    g.setColor( Constants.DOMAIN_STRUCTURE_FONT_COLOR );
                    PhylogenyDataUtil.drawString( d.getName(), xa, y1
                            + tree_panel.getMainPanel().getTreeFontSet()._fm_small.getAscent() + 6, g );
                }
                drawDomain( xa, y1, xb - xa, _rendering_height, d.getName(), g, to_pdf );
            }
        }
    }

    public void setDisplayDomainNames( final boolean display_domain_names ) {
        _display_domain_names = display_domain_names;
    }

    public void setParameter( final double e_value_threshold_exp ) {
        _e_value_threshold_exp = ( int ) e_value_threshold_exp;
    }

    public void setRenderingFactorWidth( final double rendering_factor_width ) {
        _rendering_factor_width = rendering_factor_width;
    }

    public void setRenderingHeight( final double rendering_height ) {
        _rendering_height = rendering_height;
    }

    @Override
    public StringBuffer toNHX() {
        return _domain_structure.toNHX();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        _domain_structure.toPhyloXML( writer, level, indentation );
    }

    public static void setColorMap( final Map<String, Color> domain_colors ) {
        RenderableDomainArchitecture.Domain_colors = domain_colors;
    }
}

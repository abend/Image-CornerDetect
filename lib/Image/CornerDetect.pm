=head1 NAME

Image::CornerDetect - An implementation of the Xiaochen He corner
detection algorithm.

=head1 SYNOPSIS

  use Image::CornerDetect;
  my $detector = Image::CornerDetect->new();
  my @corners = $detector->process($filename);

=head1 DESCRIPTION

Perform corner detection on an input image, returning a list of corner points.

Converted from the Matlab at
http://www.mathworks.com/matlabcentral/fileexchange/7652-a-corner-detector-based-on-global-and-local-curvature-properties

=head1 USAGE

  use Image::CornerDetect;
  my $detector = Image::CornerDetect->new();
  my @corners = $detector->process($filename);

=head1 BUGS

You tell me.

=head1 AUTHOR

    Sasha Kovar
    CPAN ID: ABEND
    sasha-cpan@arcocene.org

=head1 COPYRIGHT

This program licensed under the BSD license.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

Image::EdgeDetect

=cut

package Image::CornerDetect;

use strict;
use warnings;
use 5.010;

use File::Basename;
use File::Path;
use List::AllUtils qw(:all);
use Image::Magick;
use Data::Dumper;
use Math::Trig;
use Math::Complex;
use POSIX qw(ceil floor);

BEGIN {
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION     = '0.9';
    @ISA         = qw(Exporter);
    @EXPORT      = qw();
    @EXPORT_OK   = qw();
    %EXPORT_TAGS = ();
}

our $EPSILON = .0001;


=head2 my $detector = Image::CornerDetect->new(\%args)

Create a new corner detector, passing in parameters to override the
defaults if desired.  Arguments are:

  write_debug_image (default off)

=cut

sub new {
  my ($this, $args) = @_;
  my $class = ref($this) || $this;

  my %args = ref($args) eq 'HASH' ? %$args : ();

  my $self =
  {
   write_debug_image => $args{write_debug_image},
  };

  bless $self, $class;

  return $self;
}


=head2 $detector->process($filename)

Perform corner detection on the input image, returning a list of
corner points.  ([x1, y1], [x2, y2],...)

=cut

sub process {
  my ($self, $filename) = @_;

  # takes as input a canny edge detected image
  my $image = Image::Magick->new;
  my $status = $image->Read($filename);
  die "read: $status\n" if $status;

  my ($curve, $curve_start, $curve_end, $curve_mode, $curve_num) = extract_curve($image);

  my @cout = get_corner($curve, $curve_start, $curve_end, $curve_mode, $curve_num,
                        $image);

  if ($self->{write_debug_image}) {
    # mark corners
    for my $corner (@cout) {
      my ($x, $y) = @$corner;
      $image->Draw(fill=>'red', primitive=>'rectangle', points=>sprintf("%d,%d %d,%d", $x-2, $y-2, $x+2, $y+2));
    }

    my ($name, $path, $suffix) = fileparse($filename, qr/\.[^.]*/);
    $image->Write("$path/${name}-out$suffix");
  }

  return @cout;
}


sub get_corner {
  my ($curve, $curve_start, $curve_end, $curve_mode, $curve_num, $bw) = @_;

  my @curve = @$curve;
  my @curve_start = @$curve_start;
  my @curve_end = @$curve_end;
  my @curve_mode = @$curve_mode;

  # denotes the standard deviation of the Gaussian filter when
  # computeing curvature. The default sig is 3.
  my $sig = 3;

  # a flag to control whether add the end points of a curve as corner,
  # 1 means Yes and 0 means No. The default value is 1.
  my $Endpoint = 1;

  # denotes the minimum ratio of major axis to minor axis of an
  # ellipse, whose vertex could be detected as a corner by proposed
  # detector.  The default value is 1.5.
  my $C = 1.5;

  # denotes the maximum obtuse angle that a corner can have when it is
  # detected as a true corner, default value is 162.
  my $T_angle = 162;

  my @cout;

  my $GaussianDieOff = .0001;
  my @pw = 1..30;
  my $ssq = $sig * $sig;
  my $width = max(indexes { $_ > $GaussianDieOff } map { exp(-($_ ** 2) / (2 * $ssq)) } @pw);
  $width = 0 unless $width;
  $width += 1; # 1-based indexes

  my @t = (-$width..$width);
  my @gau = map { exp(-($_ ** 2) / (2 * $ssq)) / (2 * pi * $ssq) } @t;
  @gau = map { $_ / sum(@gau) } @gau;

  for (my $i = 0; $i < $curve_num; $i++) {
    # NOTE: matlab has the concept of x and y swapped, so what is x
    # here is y in matlab.
    my @x = map { $$_[0] } @{$curve[$i]};
    my @y = map { $$_[1] } @{$curve[$i]};
    my $W = $width;
    my $L = scalar(@{$curve[$i]});

    if ($L > $W) {
      # Calculate curvature
      my (@x1, @y1);
      if ($curve_mode[$i] eq 'loop') {
        @x1 = (@x[$L-$W..$L-1], @x, @x[0..$W-1]);
        @y1 = (@y[$L-$W..$L-1], @y, @y[0..$W-1]);
       } else {
         my @x11 = map { 2 * $x[0] } 1..$W;
         my @revx = reverse(@x[1..$W]);
         @x11 = pairwise { $a - $b } @x11, @revx;
         my @x12 = map { 2 * $x[$L-1] } 1..$W;
         @revx = reverse(@x[$L-$W-1..$L-2]);
         @x12 = pairwise { $a - $b } @x12, @revx;
         @x1 = (@x11, @x, @x12);

         my @y11 = map { 2 * $y[0] } 1..$W;
         my @revy = reverse(@y[1..$W]);
         @y11 = pairwise { $a - $b } @y11, @revy;
         my @y12 = map { 2 * $y[$L-1] } 1..$W;
         @revy = reverse(@y[$L-$W-1..$L-2]);
         @y12 = pairwise { $a - $b } @y12, @revy;
         @y1 = (@y11, @y, @y12);
       }

      my @xx = conv(\@x1, \@gau);
      @xx = @xx[$W..$L+3*$W-1];
      my @yy = conv(\@y1, \@gau);
      @yy = @yy[$W..$L+3*$W-1];

      my @Xu_in1 = @xx[2..$L+2*$W-1];
      my @Xu_in2 = @xx[0..$L+2*$W-3];
      my @Xu_in = map { ($Xu_in1[$_] - $Xu_in2[$_]) / 2 } 0..$#Xu_in1;
      my @Xu = ($xx[1]-$xx[0], @Xu_in, $xx[$L+2*$W-1]-$xx[$L+2*$W-2]);

      my @Yu_in1 = @yy[2..$L+2*$W-1];
      my @Yu_in2 = @yy[0..$L+2*$W-3];
      my @Yu_in = map { ($Yu_in1[$_] - $Yu_in2[$_]) / 2 } 0..$#Yu_in1;
      my @Yu = ($yy[1]-$yy[0], @Yu_in, $yy[$L+2*$W-1]-$yy[$L+2*$W-2]);

      my @Xuu_in1 = @Xu[2..$L+2*$W-1];
      my @Xuu_in2 = @Xu[0..$L+2*$W-3];
      my @Xuu_in = map { ($Xuu_in1[$_] - $Xuu_in2[$_]) / 2 } 0..$#Xuu_in1;
      my @Xuu = ($Xu[1]-$Xu[0], @Xuu_in, $Xu[$L+2*$W-1]-$Xu[$L+2*$W-2]);

      my @Yuu_in1 = @Yu[2..$L+2*$W-1];
      my @Yuu_in2 = @Yu[0..$L+2*$W-3];
      my @Yuu_in = map { ($Yuu_in1[$_] - $Yuu_in2[$_]) / 2 } 0..$#Yuu_in1;
      my @Yuu = ($Yu[1]-$Yu[0], @Yuu_in, $Yu[$L+2*$W-1]-$Yu[$L+2*$W-2]);

      my @K;
      for (my $i = 0; $i < scalar @Xu; $i++) {
        my $k1 = $Xu[$i] * $Yuu[$i] - $Xuu[$i] * $Yu[$i];
        my $k2 = ($Xu[$i] * $Xu[$i] + $Yu[$i] * $Yu[$i]) ** 1.5;
        push @K, abs($k1 / $k2);
      }
      @K = map { ceil($_ * 100) / 100 } @K;

      # Find curvature local maxima as corner candidates
      my @extremum;
      my $N = scalar @K;
      my $Search = 1;

      for (my $j = 0; $j < $N - 1; $j++) {
        if (($K[$j+1] - $K[$j]) * $Search > 0) {
          push @extremum, $j; # In extremum, odd points is minima and even points is maxima
          $Search = -$Search;
        }
      }

      # add a last min if we need one
      if (@extremum % 2 == 0) {
        push @extremum, $N - 1;
      }

      my @flag = map { 1 } 1..@extremum;

      # Compare with adaptive local threshold to remove round corners
      for (my $j = 0; $j < @extremum - 2; $j += 2) {
        my @tmp = @K[reverse($extremum[$j]..$extremum[$j+1])];
        my $x = min(@tmp);
        my $index1 = first_index { $_ == $x } @tmp;

        @tmp = @K[$extremum[$j+1]..$extremum[$j+2]];
        $x = min(@tmp);
        my $index2 = first_index { $_ == $x } @tmp;

        my @ROS = @K[$extremum[$j+1]-$index1..$extremum[$j+1]+$index2];

        my $K_thre = $C * mean(@ROS);
        if ($K[$extremum[$j]] < $K_thre) {
          $flag[$j] = 0;
        }
      }

      my $j = 0;
      @extremum = grep { $j++ % 2 } @extremum;
      $j = 0;
      @flag = grep { $j++ % 2 } @flag;

      @extremum = @extremum[indexes { $_ == 1 } @flag];

      # Check corner angle to remove false corners due to boundary noise and trivial details
      @flag = (0);
      our ($a, $b); # hush warnings
      my @smoothed_curve = pairwise { [$b,$a] } @xx, @yy; # NOTE: swapped x and y

      while (scalar(grep { $_ == 0 } @flag) > 0) {
        my $n = scalar @extremum;
        @flag = map { 1 } 1..@extremum;

        for (my $j = 0; $j < $n; $j++) {
          my $ang;
          if ($j == 0 && $j == $n - 1) {
            my @points = @smoothed_curve[0..min($#smoothed_curve, $L+2*$W-1)];
            $ang = curve_tangent(\@points, $extremum[$j]);
          } elsif ($j == 0) {
            my @points = @smoothed_curve[0..$extremum[$j+1]];
            $ang = curve_tangent(\@points, $extremum[$j]+1);
          } elsif ($j == $n-1) {
            my @points = @smoothed_curve[$extremum[$j-1]..min($#smoothed_curve, $L+2*$W-1)];
            my $center = $extremum[$j]-$extremum[$j-1]+1;
            $ang = curve_tangent(\@points, $center);
          } else {
            my @points = @smoothed_curve[$extremum[$j-1]..$extremum[$j+1]];
            my $center = $extremum[$j]-$extremum[$j-1]+1;
            $ang = curve_tangent(\@points, $center);
          }
          if ($ang > $T_angle && $ang < (360 - $T_angle)) {
            $flag[$j] = 0;
          }
        }

        if (scalar @extremum == 0) {
          @extremum = ();
        } else {
          @extremum = @extremum[indexes { $_ != 0 } @flag];
        }
      }

      @extremum = map { $_ - $W } @extremum;
      @extremum = grep { $_ >= 0 && $_ < $L } @extremum;
      for my $j (0..$#extremum) {
        push @cout, $curve[$i][$extremum[$j]];
      }
    }
  }

  # Add Endpoints
  if ($Endpoint) {
    for my $i (0..$#curve) {
      if (scalar(@{$curve[$i]}) > 0 && $curve_mode[$i] eq 'line') {
        # Start point compare with detected corners
        my @start_point = @{$curve_start[$i]};
        my @compare_corner = map { ($$_[0] - $start_point[0]) ** 2 +
                                     ($$_[1] - $start_point[1]) ** 2 } @cout;
        # Add end points far from detected corners
        if (min(@compare_corner) > 25) {
          push @cout, \@start_point;
        }

        # End point compare with detected corners
        my @end_point = @{$curve_end[$i]};
        @compare_corner = map { ($$_[0] - $end_point[0]) ** 2 +
                                  ($$_[1] - $end_point[1]) ** 2 } @cout;
        # Add end points far from detected corners
        if (min(@compare_corner) > 25) {
          push @cout, \@end_point;
        }
      }
    }
  }

  return @cout;
}

# $cur is an arrayref of points.
# $center is an index into cur of the middle of the curve.
sub curve_tangent {
  my ($cur, $center) = @_;
  my @cur = @$cur;

  my @direction = (0, 0);

  for my $i (1..2) {

    # once for each half of the curve
    my @curve;
    if ($i == 1) {
      @curve = reverse @cur[0..$center-1];
    } else {
      @curve = @cur[$center-1..$#cur];
    }

    my $L = scalar(@curve);

    my $tangent_direction;
    if ($L > 3) {
      my ($x1, $y1, $x2, $y2, $x3, $y3);
      # endpoints not the same?
      if (abs($curve[0][0] - $curve[$L-1][0]) > $EPSILON or
          abs($curve[0][1] - $curve[$L-1][1]) > $EPSILON)
      {
        # choose start, mid, endpoints
        my $M = ceil($L / 2); # midpoint index
        $x1 = $curve[0][0];
        $y1 = $curve[0][1];
        $x2 = $curve[$M-1][0];
        $y2 = $curve[$M-1][1];
        $x3 = $curve[$L-1][0];
        $y3 = $curve[$L-1][1];
      } else {
        # loop: choose start, 1/3, and 2/3 through
        my $M1 = ceil($L / 3);
        my $M2 = ceil(2 * $L / 3);
        $x1 = $curve[0][0];
        $y1 = $curve[0][1];
        $x2 = $curve[$M1-1][0];
        $y2 = $curve[$M1-1][1];
        $x3 = $curve[$M2-1][0];
        $y3 = $curve[$M2-1][1];
      }

      if (abs(($x1-$x2)*($y1-$y3)-($x1-$x3)*($y1-$y2))<1e-8) {
        # straight line
        $tangent_direction = theta(cplx($curve[$L-1][0] - $curve[0][0], $curve[$L-1][1] - $curve[0][1]));
      } else {
        # Fit a circle
        my $x0 = 1/2*(-$y1*$x2**2+$y3*$x2**2-$y3*$y1**2-$y3*$x1**2-$y2*$y3**2+$x3**2*$y1+$y2*$y1**2-$y2*$x3**2-$y2**2*$y1+$y2*$x1**2+$y3**2*$y1+$y2**2*$y3)/(-$y1*$x2+$y1*$x3+$y3*$x2+$x1*$y2-$x1*$y3-$x3*$y2);
        my $y0 = -1/2*($x1**2*$x2-$x1**2*$x3+$y1**2*$x2-$y1**2*$x3+$x1*$x3**2-$x1*$x2**2-$x3**2*$x2-$y3**2*$x2+$x3*$y2**2+$x1*$y3**2-$x1*$y2**2+$x3*$x2**2)/(-$y1*$x2+$y1*$x3+$y3*$x2+$x1*$y2-$x1*$y3-$x3*$y2);
        # R = ($x0-$x1)**2+($y0-$y1)**2;

        my $cp1 = cplx($x0 - $x1, $y0 - $y1);
        my $radius_direction = theta(cplx($x0 - $x1, $y0 - $y1));
        my $adjacent_direction = theta(cplx($x2 - $x1, $y2 - $y1));
        my $sin = sin($adjacent_direction - $radius_direction);
        my $sign = $sin == 0 ? 0 : $sin > 0 ? 1 : -1;
        $tangent_direction = $sign * pi / 2 + $radius_direction;
      }
    } else { # very short line
      $tangent_direction = theta(cplx($curve[$L-1][0] - $curve[0][0], $curve[$L-1][1] - $curve[0][1]));
    }
    $direction[$i-1] = $tangent_direction * 180/pi;
  }

  my $ang = abs($direction[0] - $direction[1]);

  return $ang;
}

# Function to extract curves from binary edge map, if the endpoint of a
# contour is nearly connected to another endpoint, fill the gap and continue
# the extraction. The default gap size is 1 pixels.
sub extract_curve {
  my ($bw, $gap_size) = @_;

  $gap_size ||= 1;

  my $l = $bw->Get('height');
  my $w = $bw->Get('width');

  my $bw1 = Image::Magick->new;
  $bw1->Set(size=>($w+2*$gap_size).'x'.($l+2*$gap_size));
  $bw1->ReadImage('xc:black');

  my $bw_edge = Image::Magick->new;
  $bw_edge->Set(size=>$w.'x'.$l);
  $bw_edge->ReadImage('xc:black');

  $bw1->Composite(image=>$bw,
                  compose=>'over',
                  gravity=>'center');

  # indices of nonzero values in BW1
  my @curves;
  my $FIRST =0;
  for (my @points = find_ones($bw1); @points > 0; @points = find_ones($bw1)) {
    my $point = $points[0];

    my @cur = ($point);

    set_pixel($bw1, $point, 0);

    for (my @IJ = find_ones($bw1, $point, $gap_size);
         @IJ > 0;
         @IJ = find_ones($bw1, $point, $gap_size))
    {
      my @dist;
      for my $p (@IJ) {
        my $px = ($$p[0] - $gap_size - 1) ** 2;
        my $py = ($$p[1] - $gap_size - 1) ** 2;
        my $dist = $px + $py;
        push @dist, $dist;
      }

      my $index = min_idx(@dist);

      $point = [$$point[0] + $IJ[$index][0] - $gap_size - 1,
                $$point[1] + $IJ[$index][1] - $gap_size - 1];

      push @cur, $point;

      set_pixel($bw1, $point, 0);
    }

    # Extract edge towards another direction
    $point = $points[0];

    set_pixel($bw1, $point, 0);

    for (my @IJ = find_ones($bw1, $point, $gap_size);
         @IJ > 0;
         @IJ = find_ones($bw1, $point, $gap_size))
    {
      my @dist;
      for my $p (@IJ) {
        my $px = ($$p[0] - $gap_size - 1) ** 2;
        my $py = ($$p[1] - $gap_size - 1) ** 2;
        my $dist = $px + $py;
        push @dist, $dist;
      }

      my $index = min_idx(@dist);

      $point = [$$point[0] + $IJ[$index][0] - $gap_size - 1,
                $$point[1] + $IJ[$index][1] - $gap_size - 1];

      unshift @cur, $point;

      set_pixel($bw1, $point, 0);
    }

    if (@cur > ($l + $w) / 25) {  # curve large enough? add to list.
      for my $c (@cur) {
        $$c[0] -= $gap_size;
        $$c[1] -= $gap_size;
      }
      push @curves, [@cur];
    } else {
      #say "not enough curve points (".@cur.") for ".(($l + $w) / 25);
    }
  }


  my @curve_start;
  my @curve_end;
  my @curve_mode;

  for my $curve (@curves) {
    my $curve_start = $$curve[0];
    my $curve_end = $$curve[scalar(@$curve)-1];

    my $curve_mode = 'line';
    if ((($$curve_start[0] - $$curve_end[0]) ** 2 +
        ($$curve_start[1] - $$curve_end[1]) ** 2) <= 32)
    {
      $curve_mode = 'loop';
    }

    push @curve_start, $curve_start;
    push @curve_end, $curve_end;
    push @curve_mode, $curve_mode;

    for my $p (@$curve) {
      set_pixel($bw_edge, $p, 1);
    }
  }

  #$bw_edge->write('edge.gif');

  return \@curves, \@curve_start, \@curve_end, \@curve_mode, scalar @curves;
}

sub find_ones {
  my ($img, $point, $gap) = @_;

  my ($x1, $x2, $y1, $y2, $xoff, $yoff);

  if ($point) {
    $x1 = $$point[0] - $gap;
    $x2 = $$point[0] + $gap;
    $y1 = $$point[1] - $gap;
    $y2 = $$point[1] + $gap;
    $xoff = $x1 - 1;
    $yoff = $y1 - 1;
  } else {
    # default to full image if no subrange specified
    my $w = $img->Get('width');
    $x1 = 0;
    $x2 = $w - 1;
    $xoff = 0;

    my $h = $img->Get('height');
    $y1 = 0;
    $y2 = $h - 1;
    $yoff = 0;
  }

  my @points;

  for (my $x = $x1; $x <= $x2; $x++) {
    for (my $y = $y1; $y <= $y2; $y++) {
      if (get_pixel($img, $x, $y) > .7) { # NOTE cutoff
        push @points, [$x - $xoff, $y - $yoff];
      }
    }
  }

  return @points;
}

sub get_pixel {
  my ($img, $x, $y) = @_;

  # convert from 1-based to 0-based index
  $x -= 1;
  $y -= 1;

  my @c = $img->GetPixel(x => $x, y => $y);
  return $c[0];
}

sub set_pixel {
  my ($img, $point, $val) = @_;

  my ($x, $y) = @$point;

  # convert from 1-based to 0-based index
  $x -= 1;
  $y -= 1;

  $img->SetPixel(x => $x, y => $y, color => [$val, $val, $val]);
}

sub min_idx {
  my @vals = @_;

  my $idx;
  my $cur;

  for (my $i = 0; $i < @vals; $i++) {
    my $v = $vals[$i];

    if (!$cur or $v < $cur) {
      $idx = $i;
      $cur = $v;
    }
  }

  return $idx;
}

# convolve two vectors, as http://www.mathworks.com/help/techdoc/ref/conv.html
sub conv {
  my ($a, $b) = @_;
  my @a = @$a;
  my @b = @$b;

  my @r;
  for my $x (0 .. $#a) {
    my ($col, $mult) = ($x, $a[$x]);
    for my $y (0 .. $#b) {
      $r[$col++] += $b[$y] * $mult;
    }
  }

  return @r;
}

sub mean { return sum(@_) / @_; }

1;

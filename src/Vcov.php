<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Formula\Data;
use Phit\Polyfill\Internal\Statrs;

/** Sandwich estimators from public mathematical specifications, without GPL source. */
final class Vcov
{
    private static function invalid(string $message): \Exception { return ErrorCode::exception(1001,$message); }

    public static function spec(mixed $vcov,mixed $data): array
    {
        if ($vcov===null) throw new \TypeError('Argument `$vcov` must not be null');
        $spec=['method'=>$vcov,'cluster'=>null,'label'=>'iid'];
        if ($vcov==='iid') return $spec;
        if (is_string($vcov)) {
            if (!in_array($vcov,['HC0','HC1','HC2','HC3'],true)) throw self::invalid('vcov must be iid, HC0-HC3, cluster or newey_west');
            $spec['label']=$vcov;return $spec;
        }
        if (!is_array($vcov)) throw self::invalid('vcov must be a string or an options array');
        if (array_key_exists('cluster',$vcov)) {
            if (array_diff(array_keys($vcov),['cluster','type'])!==[]) throw self::invalid('cluster options accept only cluster and type');
            $type=$vcov['type']??'HC1';
            if (!in_array($type,['HC0','HC1'],true) || (array_key_exists('type',$vcov) && $vcov['type']===null)) throw self::invalid('cluster type must be HC0 or HC1');
            $frame=Data::from($data);$input=$vcov['cluster'];
            if (is_string($input)) {
                $column=$frame->column($input);
                if ($column===null || $column['t']==='invalid') throw self::invalid('unknown or invalid cluster column');
                $values=$column['t']==='factor' ? $column['codes'] : $column['v'];
            } elseif (is_array($input)) $values=array_values($input);
            else throw self::invalid('cluster must be a column name or array of labels');
            if (count($values)!==$frame->rows) throw ErrorCode::exception(1002,'cluster labels must align with original rows');
            $ids=[];$labels=[];
            foreach ($values as $value) {
                if ($value===null || (is_float($value) && !is_finite($value))) {$labels[]=null;continue;}
                if (!is_scalar($value)) throw self::invalid('cluster labels must be scalar values or null');
                $key=is_string($value)?'s:'.$value:(is_bool($value)?'b:'.(int)$value:'n:'.(is_int($value)?(string)$value:((float)(int)$value===$value?(string)(int)$value:sprintf('%.17g',$value))));
                if (!array_key_exists($key,$ids)) $ids[$key]=count($ids);
                $labels[]=(float)$ids[$key];
            }
            $spec['cluster']=$labels;$spec['label']='cluster ('.$type.', '.(is_string($input)?$input:'supplied labels').')';
        } elseif (array_key_exists('newey_west',$vcov)) {
            if (array_diff(array_keys($vcov),['newey_west','adjust'])!==[]) throw self::invalid('Newey-West options accept only newey_west and adjust');
            if (!is_int($vcov['newey_west']) || $vcov['newey_west']<0) throw self::invalid('Newey-West lag must be a nonnegative integer');
            if (array_key_exists('adjust',$vcov) && !is_bool($vcov['adjust'])) throw self::invalid('adjust must be boolean');
            $spec['label']='Newey-West (lag '.$vcov['newey_west'].', adjust '.(($vcov['adjust']??false)?'true':'false').')';
        } else throw self::invalid('vcov options require cluster or newey_west');
        return $spec;
    }

    public static function apply(array $fit,object $design,?array $z,array $weights,array $spec): array
    {
        $method=$spec['method'];$df=$fit['df_residual'];
        if ($method!=='iid') {
            [$covariance,$df]=self::covariance($design->x,$z,$weights,$fit['residuals'],$method,$spec['cluster']===null?null:array_map(fn($v)=>(int)$v,$design->extra[array_key_last($design->extra)]));
            $se=$statistics=$pValues=$ci=[];$critical=Statrs::tCritical95($df);
            foreach ($fit['coefficients'] as $j=>$b) {
                $s=sqrt($covariance[$j][$j]);$t=fdiv($b,$s);$se[]=$s;$statistics[]=$t;$pValues[]=2*Statrs::tSf(abs($t),$df);$ci[]=[$b-$critical*$s,$b+$critical*$s];
            }
            $fit=array_replace($fit,['covariance'=>$covariance,'std_errors'=>$se,'t_statistics'=>$statistics,'p_values'=>$pValues,'conf_int'=>$ci]);
        }
        return array_replace($fit,['vcov'=>$method,'inference_df'=>$df,'vcov_label'=>$spec['label']]);
    }

    private static function dot(array $a,array $b): float { $sum=0.0;foreach($a as $i=>$v)$sum+=$v*$b[$i];return $sum; }

    public static function covariance(array $x,?array $z,array $weights,array $residuals,mixed $method,?array $groups): array
    {
        $n=count($x);$p=count($x[0]);$max=max($weights);$root=array_map(fn($w)=>sqrt($w/$max),$weights);
        $design=array_map(fn($row,$r)=>array_map(fn($v)=>$v*$r,$row),$x,$root);
        if ($z!==null) {
            $wz=array_map(fn($row,$r)=>array_map(fn($v)=>$v*$r,$row),$z,$root);$projected=array_fill(0,$n,[]);
            for($j=0;$j<$p;$j++) { $fitted=Iv::leastSquares($wz,array_column($design,$j))[2];foreach($fitted as $i=>$v)$projected[$i][]=$v; }
            $design=$projected;
        }
        $bread=Iv::leastSquares($design,array_fill(0,$n,0.0))[1];$influence=[];
        foreach($design as $i=>$row) {
            $effect=[];for($j=0;$j<$p;$j++)$effect[]=self::dot($row,array_column($bread,$j));
            $scale=$root[$i]*$residuals[$i];
            if (is_string($method)) {
                $h=self::dot($row,$effect);
                if (in_array($method,['HC2','HC3'],true) && 1-$h<=16*PHP_FLOAT_EPSILON) throw ErrorCode::exception(1005,'HC2/HC3 undefined at unit leverage');
                $scale*=match($method){'HC1'=>sqrt($n/($n-$p)),'HC2'=>1/sqrt(1-$h),'HC3'=>1/(1-$h),default=>1.0};
            }
            $influence[]=array_map(fn($v)=>$scale*$v,$effect);
        }
        $df=$n-$p;$factor=1.0;
        if (is_array($method) && array_key_exists('cluster',$method)) {
            $summed=[];foreach($groups as $i=>$g) {if(!isset($summed[$g]))$summed[$g]=array_fill(0,$p,0.0);foreach($influence[$i] as $j=>$v)$summed[$g][$j]+=$v;}
            $g=count($summed);if($g<2)throw self::invalid('cluster covariance requires at least two clusters');
            $df=$g-1;$factor=$g/($g-1)*(($method['type']??'HC1')==='HC1'?($n-1)/($n-$p):1);$influence=array_values($summed);
        }
        $cov=array_fill(0,$p,array_fill(0,$p,0.0));
        foreach($influence as $row) for($j=0;$j<$p;$j++)for($k=0;$k<$p;$k++)$cov[$j][$k]+=$row[$j]*$row[$k];
        if (is_array($method) && array_key_exists('newey_west',$method)) {
            $lag=$method['newey_west'];if($lag>=$n)throw self::invalid('Newey-West lag must be less than the retained row count');
            for($h=1;$h<=$lag;$h++){ $weight=1-$h/($lag+1);for($i=$h;$i<$n;$i++)for($j=0;$j<$p;$j++)for($k=0;$k<$p;$k++)$cov[$j][$k]+=$weight*($influence[$i][$j]*$influence[$i-$h][$k]+$influence[$i-$h][$j]*$influence[$i][$k]); }
            if($method['adjust']??false)$factor=$n/($n-$p);
        }
        foreach($cov as $j=>&$row)foreach($row as $k=>&$value){$value*=$factor;if(!is_finite($value)||($j===$k && $value<0))throw ErrorCode::exception(1005,'covariance is nonfinite or has negative variance');}
        unset($row,$value);return [$cov,$df];
    }
}

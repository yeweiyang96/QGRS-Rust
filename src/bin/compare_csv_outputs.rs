use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, PartialEq)]
struct G4Record {
    start: u32,
    end: u32,
    length: u32,
    tetrads: u32,
    y1: u32,
    y2: u32,
    y3: u32,
    gscore: u32,
    sequence: String,
}

fn parse_csv_file(path: &Path) -> Result<Vec<G4Record>, Box<dyn std::error::Error>> {
    let content = fs::read_to_string(path)?;
    let mut records = Vec::new();

    for (idx, line) in content.lines().enumerate() {
        if idx == 0 {
            // 跳过表头
            continue;
        }

        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() != 9 {
            eprintln!("⚠️  跳过格式错误的行 {}: {}", idx + 1, line);
            continue;
        }

        records.push(G4Record {
            start: parts[0].parse()?,
            end: parts[1].parse()?,
            length: parts[2].parse()?,
            tetrads: parts[3].parse()?,
            y1: parts[4].parse()?,
            y2: parts[5].parse()?,
            y3: parts[6].parse()?,
            gscore: parts[7].parse()?,
            sequence: parts[8].to_string(),
        });
    }

    Ok(records)
}

fn compare_records(mmap_records: &[G4Record], stream_records: &[G4Record]) -> (usize, Vec<String>) {
    let mut mismatches = 0;
    let mut details = Vec::new();

    if mmap_records.len() != stream_records.len() {
        details.push(format!(
            "  ⚠️  记录数量不一致: mmap={}, stream={}",
            mmap_records.len(),
            stream_records.len()
        ));
        mismatches += 1;
    }

    let min_len = mmap_records.len().min(stream_records.len());
    for i in 0..min_len {
        let mmap = &mmap_records[i];
        let stream = &stream_records[i];

        if mmap != stream {
            mismatches += 1;
            if mismatches <= 10 {
                details.push(format!("  ⚠️  第 {} 条记录不匹配:", i + 1));

                if mmap.start != stream.start || mmap.end != stream.end {
                    details.push(format!(
                        "      位置: mmap={}..{}, stream={}..{}",
                        mmap.start, mmap.end, stream.start, stream.end
                    ));
                }
                if mmap.sequence != stream.sequence {
                    details.push(format!(
                        "      序列: mmap={}, stream={}",
                        mmap.sequence, stream.sequence
                    ));
                }
                if mmap.tetrads != stream.tetrads {
                    details.push(format!(
                        "      四联体: mmap={}, stream={}",
                        mmap.tetrads, stream.tetrads
                    ));
                }
                if mmap.gscore != stream.gscore {
                    details.push(format!(
                        "      gscore: mmap={}, stream={}",
                        mmap.gscore, stream.gscore
                    ));
                }
                if mmap.y1 != stream.y1 || mmap.y2 != stream.y2 || mmap.y3 != stream.y3 {
                    details.push(format!(
                        "      间隔: mmap=({},{},{}), stream=({},{},{})",
                        mmap.y1, mmap.y2, mmap.y3, stream.y1, stream.y2, stream.y3
                    ));
                }
            }
        }
    }

    if mismatches > 10 {
        details.push(format!("  ... (省略其余 {} 处差异)", mismatches - 10));
    }

    (mismatches, details)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let (mmap_dir, stream_dir) = if args.len() >= 3 {
        (PathBuf::from(&args[1]), PathBuf::from(&args[2]))
    } else {
        (
            PathBuf::from("output/dme/mmap"),
            PathBuf::from("output/dme/stream"),
        )
    };

    if !mmap_dir.exists() {
        eprintln!("❌ mmap 目录不存在: {:?}", mmap_dir);
        std::process::exit(1);
    }

    if !stream_dir.exists() {
        eprintln!("❌ stream 目录不存在: {:?}", stream_dir);
        std::process::exit(1);
    }

    println!("════════════════════════════════════════════════════════");
    println!("🔍 比较 Mmap 和 Stream 输出文件");
    println!("════════════════════════════════════════════════════════");
    println!("Mmap 目录:   {}", mmap_dir.display());
    println!("Stream 目录: {}", stream_dir.display());
    println!("════════════════════════════════════════════════════════\n");

    // 获取 mmap 目录中的所有 CSV 文件
    let mmap_files: HashSet<String> = fs::read_dir(&mmap_dir)
        .expect("无法读取 mmap 目录")
        .filter_map(|entry| {
            let entry = entry.ok()?;
            let path = entry.path();
            if path.extension()? == "csv" {
                Some(entry.file_name().to_string_lossy().to_string())
            } else {
                None
            }
        })
        .collect();

    // 获取 stream 目录中的所有 CSV 文件
    let stream_files: HashSet<String> = fs::read_dir(&stream_dir)
        .expect("无法读取 stream 目录")
        .filter_map(|entry| {
            let entry = entry.ok()?;
            let path = entry.path();
            if path.extension()? == "csv" {
                Some(entry.file_name().to_string_lossy().to_string())
            } else {
                None
            }
        })
        .collect();

    // 找出共同的文件
    let mut common_files: Vec<String> = mmap_files.intersection(&stream_files).cloned().collect();
    common_files.sort();

    // 找出只在一个目录中的文件
    let only_mmap: Vec<String> = mmap_files.difference(&stream_files).cloned().collect();
    let only_stream: Vec<String> = stream_files.difference(&mmap_files).cloned().collect();

    if !only_mmap.is_empty() {
        println!("⚠️  只在 mmap 目录中的文件:");
        for f in &only_mmap {
            println!("    {}", f);
        }
        println!();
    }

    if !only_stream.is_empty() {
        println!("⚠️  只在 stream 目录中的文件:");
        for f in &only_stream {
            println!("    {}", f);
        }
        println!();
    }

    println!("📁 找到 {} 个共同文件\n", common_files.len());

    // 比较每个文件
    let mut total_mismatches = 0;
    let mut file_results = Vec::new();

    for file_name in &common_files {
        let mmap_path = mmap_dir.join(file_name);
        let stream_path = stream_dir.join(file_name);

        print!("🔍 比较 {}... ", file_name);

        let mmap_records = match parse_csv_file(&mmap_path) {
            Ok(r) => r,
            Err(e) => {
                println!("❌");
                eprintln!("  读取 mmap 文件失败: {}", e);
                total_mismatches += 1;
                continue;
            }
        };

        let stream_records = match parse_csv_file(&stream_path) {
            Ok(r) => r,
            Err(e) => {
                println!("❌");
                eprintln!("  读取 stream 文件失败: {}", e);
                total_mismatches += 1;
                continue;
            }
        };

        let (mismatches, details) = compare_records(&mmap_records, &stream_records);

        if mismatches == 0 {
            println!("✅ ({} 条记录)", mmap_records.len());
            file_results.push((file_name.clone(), true, mmap_records.len(), details));
        } else {
            println!("❌ ({} 处差异)", mismatches);
            file_results.push((file_name.clone(), false, mismatches, details));
            total_mismatches += 1;
        }
    }

    // 显示详细差异
    println!("\n════════════════════════════════════════════════════════");
    println!("📊 详细结果:");
    println!("════════════════════════════════════════════════════════\n");

    for (file_name, is_match, count, details) in &file_results {
        if *is_match {
            println!("✅ {}: {} 条记录完全匹配", file_name, count);
        } else {
            println!("❌ {}: {} 处差异", file_name, count);
            for detail in details {
                println!("{}", detail);
            }
            println!();
        }
    }

    // 总结
    println!("════════════════════════════════════════════════════════");
    println!("📈 总结:");
    println!("════════════════════════════════════════════════════════");
    println!("共同文件数: {}", common_files.len());
    println!("完全匹配: {}", common_files.len() - total_mismatches);
    println!("有差异: {}", total_mismatches);

    if !only_mmap.is_empty() {
        println!("只在 mmap: {}", only_mmap.len());
    }
    if !only_stream.is_empty() {
        println!("只在 stream: {}", only_stream.len());
    }

    if total_mismatches == 0 && only_mmap.is_empty() && only_stream.is_empty() {
        println!("\n✅ 所有文件完全一致!");
    } else {
        println!("\n❌ 发现差异!");
        std::process::exit(1);
    }
}
